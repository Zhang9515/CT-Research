#include "FDK.h"
// 2018/04/20 apply GPU acceleration
#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;

__global__ void GFunction(float *dev_G, const float PInt, const int Pstart, const int LP)
{

	const unsigned int Pindex = Pstart + threadIdx.x;

	// S_L filter 
	if (Pindex == 0)
		dev_G[LP - 1] = 1.0 / pow2(PI*PInt);
	else
	{
		dev_G[(LP - 1) - Pindex] = -1.0 / pow2(PInt * PI) / (4 * pow2(Pindex) - 1);
		dev_G[(LP - 1) + Pindex] = -1.0 / pow2(PInt * PI) / (4 * pow2(Pindex) - 1);
	}
}

__global__ void PreWeightFilter(float *dev_Rcov, float *dev_R, const float *dev_Pdomain, const float *dev_Xigamadomain,
	const double Distance, const float *dev_G, const float PInt, const int LP,
	const int LXigama, const int Pstart, const int Xigamastart, const int Betaindex, const int gstart,
	const int gend)
{
	const unsigned int Pindex = threadIdx.x + Pstart;
	const unsigned int Xigamaindex = blockIdx.x + Xigamastart;
	const unsigned long base_id = Betaindex * (LXigama * LP) + Xigamaindex * LP;
	const unsigned long thread_id = base_id + Pindex;

	float P = dev_Pdomain[Pindex];
	float Xigama = dev_Xigamadomain[Xigamaindex];

	double Proportion = Distance / sqrt(pow2(Distance) + pow2(P) + pow2(Xigama));

	dev_R[thread_id] = dev_R[thread_id] * Proportion;    // directly cover the input

	double Rcovsum = 0;
	// convolution

	for (int g = gstart; g < gend; g++)
	{
		//convolution
		Rcovsum += dev_R[base_id + g] * dev_G[(LP - 1) + Pindex - g];
	}

	dev_Rcov[thread_id] += PInt * Rcovsum;
}

__global__ void BackProjection(const float *dev_Rcov, float *dev_Display, const double *dev_Size,
	const int t_length, const int s_length, const int z_length, const float Beta,
	const double Distance,
	const float *dev_Pdomain, const float *dev_Xigamadomain, const float PInt, const float XigamaInt,
	const float BetaScanInt, const float minP, const float maxP, const float minXigama, const float maxXigama,
	const int betaIndex, const int LP, const int LXigama, const int T_start, const int S_start)
{
	// initialize
	const double Resolution_t = 1.0 * dev_Size[0] / t_length;
	const double Resolution_s = 1.0 * dev_Size[1] / s_length;
	const double Resolution_z = 1.0 * dev_Size[2] / z_length;
	// rotation center
	double Center_t = dev_Size[0] / 2;
	double Center_s = dev_Size[1] / 2;
	double Center_z = dev_Size[2] / 2;
	// index 
	const unsigned int Tindex = T_start + threadIdx.x;
	const unsigned int Sindex = S_start + blockIdx.x;

	double image_t, image_s, image_z, dect_t, dect_s, dect_z;
	unsigned int XigamaN1index, XigamaN2index, PN1index, PN2index;
	double P_domain1, P_domain2, Xigama_domain1, Xigama_domain2;
	double Xig1, Xig2, P1, P2;
	double Display_pBeta = 0;
	double LengthRatio;
	double Xigama, P;
	double Weight;
	unsigned long thread_id;

	for (int Zindex = 0; Zindex < z_length; Zindex++)
	{
		// this is a little different from code on MATLAB
		image_t = (Tindex + 0.5) * Resolution_t - Center_t;
		image_s = (Sindex + 0.5) * Resolution_s - Center_s;
		image_z = (Zindex + 0.5) * Resolution_z - Center_z;     // image pixel in ground coordinate

																// rotate in ground coordinate
		dect_t = image_t * cos(Beta) + image_s * sin(Beta);
		dect_s = -image_t * sin(Beta) + image_s * cos(Beta);
		dect_z = image_z;

		// define the projection position on the detector
		LengthRatio = Distance / (Distance - dect_s);
		Xigama = dect_z * LengthRatio;
		P = dect_t * LengthRatio;

		if ((P >= minP) && (P < maxP) && (Xigama >= minXigama) && (Xigama < maxXigama))
		{
			XigamaN1index = floor(fabs(Xigama - dev_Xigamadomain[0]) / XigamaInt);
			XigamaN2index = XigamaN1index + 1;
			PN1index = floor(fabs(P - dev_Pdomain[0]) / PInt);
			PN2index = PN1index + 1;

			P_domain1 = dev_Pdomain[PN1index]; P_domain2 = dev_Pdomain[PN2index];
			Xigama_domain1 = dev_Xigamadomain[XigamaN1index]; Xigama_domain2 = dev_Xigamadomain[XigamaN2index];

			//bilinear interpolation
			Xig1 = fabs(Xigama - Xigama_domain1); Xig2 = fabs(Xigama_domain2 - Xigama);
			P1 = fabs(P - P_domain1); P2 = fabs(P_domain2 - P);

			Weight = sqrt(1 + 0 * pow2(dect_z) / (pow2(Distance + dect_s) + pow2(dect_t)));
			//Display_pBeta = 1;
			Display_pBeta = (Xig2 * P2 * dev_Rcov[betaIndex * LP * LXigama + XigamaN1index * LP + PN1index]
				+ Xig1 * P2 * dev_Rcov[betaIndex * LP * LXigama + XigamaN2index * LP + PN1index]
				+ Xig2 * P1 * dev_Rcov[betaIndex * LP * LXigama + XigamaN1index * LP + PN2index]
				+ Xig1 * P1 * dev_Rcov[betaIndex * LP * LXigama + XigamaN2index * LP + PN2index])
				/ (PInt * XigamaInt) * pow2(LengthRatio) * BetaScanInt * Weight;
		}
		thread_id = Zindex * (t_length * s_length) + Sindex * t_length + Tindex;
		dev_Display[thread_id] += Display_pBeta;
	}

}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t FDKpro(float *Display, const float *R, const float *Xigamadomain, const float *Pdomain,
	const float *BetaScanRange, const double Distance, const int LBeta, const int LP, const int LXigama,
	const double *Size, const int t_length, const int s_length, const int z_length)
{
	mexPrintf("Hello FDK!\n");
	float *dev_R = 0, *dev_Rcov = 0;
	float *dev_G = 0;
	float *dev_BetaScanRange = 0, *dev_Pdomain = 0, *dev_Xigamadomain = 0;
	double *dev_Size = 0;

	float *dev_Display = 0;
	float PInt = fabs(Pdomain[1] - Pdomain[0]);
	float XigamaInt = fabs(Xigamadomain[1] - Xigamadomain[0]);
	float BetaScanInt = fabs(BetaScanRange[1] - BetaScanRange[0]);

	float maxP = MAX(Pdomain[0], Pdomain[LP - 1]);
	float minP = MIN(Pdomain[0], Pdomain[LP - 1]);
	float maxXigama = MAX(Xigamadomain[0], Xigamadomain[LXigama - 1]);
	float minXigama = MIN(Xigamadomain[0], Xigamadomain[LXigama - 1]);
	//mexPrintf("%lf %lf %lf %lf \n", maxGama, minGama, maxXigama, minXigama);

	const long LDisplay = t_length * s_length * z_length;
	const long LR = LP * LXigama * LBeta;
	const int LFilter = 2 * LP - 1;

	short thread_cubic_x = MIN(threadX, LP);
	short block_cubic_x = MIN(blockX, LXigama);
	short thread_cubic_Bp_x = MIN(threadX, t_length);
	short block_cubic_Bp_x = MIN(blockX, s_length);

	const dim3 thread_cubic(thread_cubic_x, 1, 1);
	const dim3 block_cubic(block_cubic_x, 1, 1);
	const dim3 thread_cubic_Bp(thread_cubic_Bp_x, 1, 1);
	const dim3 block_cubic_Bp(block_cubic_Bp_x, 1, 1);

	dim3 thread_cubic_residual(1, 1, 1);  // initial
	dim3 block_cubic_residual(1, 1, 1);  // initial
	dim3 thread_cubic_Bp_residual(1, 1, 1);  // initial
	dim3 block_cubic_Bp_residual(1, 1, 1);  // initial

	short LPResidual = LP % threadX;
	short LXigamaResidual = LXigama % blockX;
	short PTime = LP / threadX;
	short XigamaTime = LXigama / blockX;
	short Pstart = 0;
	short Xigamastart = 0;
	short TlengthResidual = t_length % threadX;
	short SlengthResidual = s_length % blockX;
	short T_Time = t_length / threadX;
	short S_Time = s_length / blockX;
	short T_start = 0;
	short S_start = 0;
	short gstart = 0, gend = 0;
	short gtime = 1;
	//mexPrintf("PTime: %d XigamaTime: %d\n", PTime, XigamaTime);

	if (LPResidual != 0)
	{
		thread_cubic_residual.x = LPResidual;
	}
	if (LXigamaResidual != 0)
	{
		block_cubic_residual.x = LXigamaResidual;
	}
	if (TlengthResidual != 0)
	{
		thread_cubic_Bp_residual.x = TlengthResidual;
	}
	if (SlengthResidual != 0)
	{
		block_cubic_Bp_residual.x = SlengthResidual;
	}

	mexPrintf("thread_cubic_x: %d block_cubic_x: %d LPResidual: %d LXigamaResidual: %d\n",
		thread_cubic_x, block_cubic_x, LPResidual, LXigamaResidual);

	cudaError_t cudaStatus;

	///////////////////////////////////////////////////////////////////////////////////////////////	
	mexPrintf("start cuda\n");

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
		mexPrintf("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	mexPrintf("call for space in GPU\n");

	// Allocate GPU buffers for three vectors (two input, one output).

	cudaStatus = cudaMalloc((void**)&dev_R, LR * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_R cudaMalloc failed!\n");
		mexPrintf("dev_R cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error1;
	}

	cudaStatus = cudaMalloc((void**)&dev_Pdomain, LP * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pdomain cudaMalloc failed!\n");
		mexPrintf("dev_Pdomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Xigamadomain, LXigama * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Xigamadomain cudaMalloc failed!\n");
		mexPrintf("dev_Xigamadomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_BetaScanRange, LBeta * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_BetaScanRange cudaMalloc failed!\n");
		mexPrintf("dev_BetaScanRange cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	//mexPrintf("copy data in CPU to GPU\n");

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_R, R, LR * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy R failed!\n");
		mexPrintf("cudaMemcpy R failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error1;
	}

	cudaStatus = cudaMemcpy(dev_Pdomain, Pdomain, LP * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Gamadomain failed!\n");
		mexPrintf("cudaMemcpy Gamadomain failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Xigamadomain, Xigamadomain, LXigama * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Xigamadomain failed!\n");
		mexPrintf("cudaMemcpy Xigamadomain failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_BetaScanRange, BetaScanRange, LBeta * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy BetaScanRange failed!\n");
		mexPrintf("cudaMemcpy BetaScanRange failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	mexPrintf("start parallel computation\n");

	mexPrintf("Preweight and Filtering\n");
	// Preweight and filter 
	// to be adapted to GPU, we limited the number of threads in each grid (threadX * blockX * LBeta)
	// generate filter
	mexPrintf("Generate filter\n");
	cudaStatus = cudaMalloc((void**)&dev_G, LFilter * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_G cudaMalloc failed!\n");
		mexPrintf("dev_G cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error1;
	}
	cudaMemset(dev_G, 0, sizeof(float));

	for (int numP = 0; numP < PTime; numP++)
	{
		Pstart = numP * threadX;
		GFunction << <1, thread_cubic_x >> > (dev_G, PInt, Pstart, LP);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "generate filter failed: %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("generate filter launch failed %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("Error happens at numP: %d\n", numP);
			goto Error;
		}
	}
	if (LPResidual != 0)
	{
		Pstart = LP - LPResidual;
		GFunction << <1, thread_cubic_residual.x >> > (dev_G, PInt, Pstart, LP);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "generate filter failed: %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("generate filter launch failed %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("Error happens at LXigamaResidual: %d\n", LXigamaResidual);
			goto Error;
		}
	}

	// store the result of convolution
	cudaStatus = cudaMalloc((void**)&dev_Rcov, LR * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Rcov cudaMalloc failed!\n");
		mexPrintf("dev_Rcov cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	cudaMemset(dev_Rcov, 0, sizeof(float));

	// convolution with limited scale
	gtime = ceil(1.0 * LP / Filterlengthlimit);
	mexPrintf("gtime: %d\n", gtime);

	for (int numB = 0; numB < LBeta; numB++)
	{
		//Beta = BetaScanRange[numB];
		for (int numP = 0; numP < PTime; numP++)
		{
			for (int numX = 0; numX < XigamaTime; numX++)
			{
				Pstart = numP * threadX;
				Xigamastart = numX * blockX;
				//mexPrintf("%d %d\n", Pstart, Xigamastart);
				for (int numg = 0; numg < gtime; numg++)
				{
					gstart = int(numg * Filterlengthlimit);
					if (numg != (gtime - 1))
					{
						gend = int((numg + 1)* Filterlengthlimit);
					}
					else
					{
						gend = LP;
					}
					/*mexPrintf("gstart: %d\n", gstart);
					mexPrintf("gend: %d\n", gend);*/
					PreWeightFilter << <block_cubic, thread_cubic >> > (dev_Rcov, dev_R, dev_Pdomain,
						dev_Xigamadomain, Distance, dev_G, PInt, LP, LXigama,
						Pstart, Xigamastart, numB, gstart, gend);
					// Check for any errors launching the kernel
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess) {
						fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Error happens at numB: %d Pstart: %d Xigamastart: %d numg: %d\n",
							numB, Pstart, Xigamastart, numg);
						goto Error;
					}
				}
			}
		}

		if (LPResidual != 0)
		{
			Pstart = LP - LPResidual;
			if (LXigamaResidual != 0)
			{
				Xigamastart = LXigama - LXigamaResidual;
				//("%d %d\n", Pstart, Xigamastart);
				for (int numg = 0; numg < gtime; numg++)
				{
					gstart = int(numg * Filterlengthlimit);
					if (numg != (gtime - 1))
					{
						gend = int((numg + 1)* Filterlengthlimit);
					}
					else
					{
						gend = LP;
					}
					/*mexPrintf("gstart: %d\n", gstart);
					mexPrintf("gend: %d\n", gend);*/
					PreWeightFilter << <block_cubic_residual, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_Pdomain,
						dev_Xigamadomain, Distance, dev_G, PInt, LP, LXigama,
						Pstart, Xigamastart, numB, gstart, gend);
					// Check for any errors launching the kernel
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess) {
						fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Error happens at numB: %d Pstart: %d Xigamastart: %d numg: %d\n",
							numB, Pstart, Xigamastart, numg);
						goto Error;
					}
				}
			}

			for (int numX = 0; numX < XigamaTime; numX++)
			{
				Xigamastart = numX * blockX;
				//("%d %d\n", Pstart, Xigamastart);
				for (int numg = 0; numg < gtime; numg++)
				{
					gstart = int(numg * Filterlengthlimit);
					if (numg != (gtime - 1))
					{
						gend = int((numg + 1)* Filterlengthlimit);
					}
					else
					{
						gend = LP;
					}
					/*mexPrintf("gstart: %d\n", gstart);
					mexPrintf("gend: %d\n", gend);*/
					PreWeightFilter << <block_cubic, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_Pdomain,
						dev_Xigamadomain, Distance, dev_G, PInt, LP, LXigama,
						Pstart, Xigamastart, numB, gstart, gend);
					// Check for any errors launching the kernel
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess) {
						fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Error happens at numB: %d Pstart: %d Xigamastart: %d numg: %d\n",
							numB, Pstart, Xigamastart, numg);
						goto Error;
					}
				}
			}
		}
		if (LXigamaResidual != 0)
		{
			Xigamastart = LXigama - LXigamaResidual;
			for (int numP = 0; numP < PTime; numP++)
			{
				Pstart = numP * threadX;
				//mexPrintf("%d %d\n", Pstart, Xigamastart);
				for (int numg = 0; numg < gtime; numg++)
				{
					gstart = int(numg * Filterlengthlimit);
					if (numg != (gtime - 1))
					{
						gend = int((numg + 1)* Filterlengthlimit);
					}
					else
					{
						gend = LP;
					}
					/*mexPrintf("gstart: %d\n", gstart);
					mexPrintf("gend: %d\n", gend);*/
					PreWeightFilter << <block_cubic_residual, thread_cubic >> > (dev_Rcov, dev_R, dev_Pdomain,
						dev_Xigamadomain, Distance, dev_G, PInt, LP, LXigama,
						Pstart, Xigamastart, numB, gstart, gend);
					// Check for any errors launching the kernel
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess) {
						fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
						mexPrintf("Error happens at numB: %d Pstart: %d Xigamastart: %d numg: %d\n",
							numB, Pstart, Xigamastart, numg);
						goto Error;
					}
				}
			}
		}
	}

Error1:
	cudaFree(dev_R);
	cudaFree(dev_G);

	mexPrintf("backprojection\n");

	cudaStatus = cudaMalloc((void**)&dev_Size, 3 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Size cudaMalloc failed!\n");
		mexPrintf("dev_Size cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Size, Size, 3 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Size failed!\n");
		mexPrintf("cudaMemcpy Size failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Display, LDisplay * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Display cudaMalloc failed!\n");
		mexPrintf("dev_Display cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	cudaMemset(dev_Display, 0, sizeof(float));

	//Backprojection  
	for (int betaIndex = 0; betaIndex < LBeta; betaIndex++)
	{
		for (int numT = 0; numT < T_Time; numT++)
		{
			for (int numS = 0; numS < S_Time; numS++)
			{
				T_start = numT * threadX;
				S_start = numS * blockX;
				BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length, z_length,
					BetaScanRange[betaIndex], Distance, dev_Pdomain, dev_Xigamadomain, PInt, XigamaInt, BetaScanInt, minP, maxP,
					minXigama, maxXigama, betaIndex, LP, LXigama, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at betaIndex: %d\n", betaIndex);
					goto Error;
				}
			}
		}
		if (TlengthResidual != 0)
		{
			T_start = t_length - TlengthResidual;
			for (int numS = 0; numS < S_Time; numS++)
			{
				S_start = numS * blockX;
				BackProjection << <block_cubic_Bp, thread_cubic_Bp_residual >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length, z_length,
					BetaScanRange[betaIndex], Distance, dev_Pdomain, dev_Xigamadomain, PInt, XigamaInt, BetaScanInt, minP, maxP,
					minXigama, maxXigama, betaIndex, LP, LXigama, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at betaIndex: %d\n", betaIndex);
					goto Error;
				}
			}
			if (SlengthResidual != 0)
			{
				S_start = s_length - SlengthResidual;
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp_residual >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length, z_length,
					BetaScanRange[betaIndex], Distance, dev_Pdomain, dev_Xigamadomain, PInt, XigamaInt, BetaScanInt, minP, maxP,
					minXigama, maxXigama, betaIndex, LP, LXigama, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at betaIndex: %d\n", betaIndex);
					goto Error;
				}
			}
		}
		if (SlengthResidual != 0)
		{
			S_start = s_length - SlengthResidual;
			for (int numT = 0; numT < T_Time; numT++)
			{
				T_start = numT * threadX;
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length, z_length,
					BetaScanRange[betaIndex], Distance, dev_Pdomain, dev_Xigamadomain, PInt, XigamaInt, BetaScanInt, minP, maxP,
					minXigama, maxXigama, betaIndex, LP, LXigama, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at betaIndex: %d\n", betaIndex);
					goto Error;
				}
			}
		}		
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		mexPrintf("cudaDeviceSynchronize returned error code %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	//Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(Display, dev_Display, LDisplay * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!\n");
		mexPrintf("cudaMemcpy dev_Display failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	cudaFree(dev_BetaScanRange);
	cudaFree(dev_Pdomain);
	cudaFree(dev_Xigamadomain);
	cudaFree(dev_Rcov);
	cudaFree(dev_Display);
	cudaFree(dev_Size);

	mexPrintf("Exit FDK\n");
	return cudaStatus;
}