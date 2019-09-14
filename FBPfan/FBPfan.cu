#include "FBPfan.h"
// 2019/09/14 apply GPU acceleration
#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;

__global__ void GFunction(float *dev_G, const float PInt, const int pstart, const int LP)
{

	const unsigned int Pindex = pstart + threadIdx.x;

	// S_L filter 
	if (Pindex == 0)
		dev_G[LP - 1] = 1.0 / pow2(PI*PInt);
	else
	{
		dev_G[(LP - 1) - Pindex] = -1.0 / pow2(PInt * PI) / (4 * pow2(Pindex) - 1);
		dev_G[(LP - 1) + Pindex] = -1.0 / pow2(PInt * PI) / (4 * pow2(Pindex) - 1);
	}
}

__global__ void Filter(float *dev_Rcov, float *dev_R, const float *dev_G, const float PInt, 
	const int LP, const int pstart, const int Betastart, const int gstart, const int gend)
{
	const unsigned int Pindex = threadIdx.x + pstart;
	const unsigned int BetaIndex = blockIdx.x + Betastart;
	const unsigned long base_id = BetaIndex * LP;
	const unsigned long thread_id = base_id + Pindex;

	double Rcovsum = 0;
	// convolution

	for (int g = gstart; g < gend; g++)
	{
		//convolution
		Rcovsum += dev_R[base_id + g] * dev_G[(LP - 1) + Pindex - g];
	}

	dev_Rcov[thread_id] += PInt * Rcovsum;
}

__global__ void BackProjection(const float *dev_Rcov, float *dev_Display, const double *dev_Size, const int T_length, const int S_length,
	const float Beta, const float *dev_Pdomain, const float PInt, const float BetaScanInt, const float minP, 
	const float maxP, const int BetaIndex, const int LP, const int LBeta, const int T_start, const int S_start, const double Rscan)
{
	// here makes clear an important note: t in lower case means the location on the detector, while T in captial case means a dimension of the axis.
	// only take effect in this function
	// initialize
	const double Resolution_t = 1.0 * dev_Size[0] / T_length;
	const double Resolution_s = 1.0 * dev_Size[1] / S_length;
	// rotation center
	double Center_t = dev_Size[0] / 2;
	double Center_s = dev_Size[1] / 2;
	// index 
	const unsigned int Tindex = T_start + threadIdx.x;
	const unsigned int Sindex = S_start + blockIdx.x;

	double image_t, image_s;
	unsigned int P1index, P2index;
	double P_domain1, P_domain2;
	double P1, P2;
	double Display_pBeta = 0;
	double P, S, U, P_domain;
	unsigned long thread_id;

	// this is a little different from code on MATLAB
	image_t = (Tindex + 0.5) * Resolution_t - Center_t;    // center of the image pixel in ground coordinate
	image_s = (Sindex + 0.5) * Resolution_s - Center_s;

	P = image_t * cos(Beta) + image_s * sin(Beta);   // rotate in ground coordinate
	S = -image_t * sin(Beta) + image_s * cos(Beta);
	U = (Rscan - S) / Rscan;       // proportion of length

	P_domain = P / U;
	
	if ((P_domain >= minP) && (P_domain < maxP))
	{
		P1index = floor(fabs(P_domain - dev_Pdomain[0]) / PInt);
		P2index = P1index + 1;

		P_domain1 = dev_Pdomain[P1index]; P_domain2 = dev_Pdomain[P2index];
	
		//bilinear interpolation
		P1 = fabs(P_domain - P_domain1); P2 = fabs(P_domain2 - P_domain);

		//Display_pBeta = 1;
		Display_pBeta = ( P2 * dev_Rcov[BetaIndex * LP + P1index] + P1 * dev_Rcov[BetaIndex * LP + P2index] )
			/ (pow2(U) * PInt) * BetaScanInt;
	}
	thread_id = Sindex * T_length + Tindex;
	dev_Display[thread_id] += Display_pBeta;

}

// HeLTer function for using CUDA to add vectors in parallel.
cudaError_t FBPfan(float *Display, const float *R, const float *Pdomain, const float *BetaScanRange, const int LBeta,
	const int LP, const double *Size, const int t_length, const int s_length, const double Rscan)
{
	mexPrintf("Hello FBPpara!\n");
	float *dev_R = 0, *dev_Rcov = 0;
	float *dev_G = 0;
	float *dev_BetaScanRange = 0, *dev_Pdomain = 0;
	double *dev_Size = 0;

	float *dev_Display = 0;
	float PInt = fabs(Pdomain[1] - Pdomain[0]);
	float BetaScanInt = fabs(BetaScanRange[1] - BetaScanRange[0]);

	float maxP = MAX(Pdomain[0], Pdomain[LP - 1]);
	float minP = MIN(Pdomain[0], Pdomain[LP - 1]);
	//mexPrintf("%lf %lf %lf %lf \n", maxGama, minGama, maxXigama, minXigama);

	const long LDisplay = t_length * s_length;
	const long LR = LP * LBeta;
	const int LFilter = 2 * LP - 1;

	short thread_cubic_x = MIN(threadX, LP);
	short block_cubic_x = MIN(blockX, LBeta);
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
	short LBetaResidual = LBeta % blockX;
	short pTime = LP / threadX;
	short BetaTime = LBeta / blockX;
	short pstart = 0;
	short Betastart = 0;
	short TlengthResidual = t_length % threadX;
	short SlengthResidual = s_length % blockX;
	short T_Time = t_length / threadX;
	short S_Time = s_length / blockX;
	short T_start = 0;
	short S_start = 0;
	short gstart = 0, gend = 0;
	short gtime = 1;
	//mexPrintf("pTime: %d BetaTime: %d\n", pTime, BetaTime);

	if (LPResidual != 0)
	{
		thread_cubic_residual.x = LPResidual;
	}
	if (LBetaResidual != 0)
	{
		block_cubic_residual.x = LBetaResidual;
	}
	if (TlengthResidual != 0)
	{
		thread_cubic_Bp_residual.x = TlengthResidual;
	}
	if (SlengthResidual != 0)
	{
		block_cubic_Bp_residual.x = SlengthResidual;
	}

	mexPrintf("thread_cubic_x: %d block_cubic_x: %d LPResidual: %d LBetaResidual: %d\n",
		thread_cubic_x, block_cubic_x, LPResidual, LBetaResidual);
	mexPrintf("thread_cubic_Bp_x: %d block_cubic_Bp_x: %d TlengthResidual: %d SlengthResidual: %d\n",
		thread_cubic_Bp_x, block_cubic_Bp_x, TlengthResidual, SlengthResidual);

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
		fprintf(stderr, "cudaMemcpy Pdomain failed!\n");
		mexPrintf("cudaMemcpy Pdomain failed! %s\n", cudaGetErrorString(cudaStatus));
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

	for (int numP = 0; numP < pTime; numP++)
	{
		pstart = numP * threadX;
		GFunction << <1, thread_cubic_x >> > (dev_G, PInt, pstart, LP);
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
		pstart = LP - LPResidual;
		GFunction << <1, thread_cubic_residual.x >> > (dev_G, PInt, pstart, LP);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "generate filter failed: %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("generate filter launch failed %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("Error happens at LPResidual: %d\n", LPResidual);
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

	for (int numP = 0; numP < pTime; numP++)
	{
		for (int numBeta = 0; numBeta < BetaTime; numBeta++)
		{
			pstart = numP * threadX;
			Betastart = numBeta * blockX;
			//mexPrintf("%d %d\n", pstart, Betastart);
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
				Filter << <block_cubic, thread_cubic >> > (dev_Rcov, dev_R, dev_G, PInt, LP, pstart, Betastart, gstart, gend);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at pstart: %d Betastart: %d numg: %d\n",
						pstart, Betastart, numg);
					goto Error;
				}
			}
		}
	}

	if (LPResidual != 0)
	{
		pstart = LP - LPResidual;
		if (LBetaResidual != 0)
		{
			Betastart = LBeta - LBetaResidual;
			//("%d %d\n", pstart, Betastart);
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
				Filter << <block_cubic_residual, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_G, PInt, LP, pstart, Betastart, gstart, gend);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at pstart: %d Betastart: %d numg: %d\n",
						pstart, Betastart, numg);
					goto Error;
				}
			}
		}

		for (int numBeta = 0; numBeta < BetaTime; numBeta++)
		{
			Betastart = numBeta * blockX;
			//("%d %d\n", pstart, Betastart);
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
				Filter << <block_cubic, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_G, PInt, LP, pstart, Betastart, gstart, gend);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at pstart: %d Betastart: %d numg: %d\n",
						pstart, Betastart, numg);
					goto Error;
				}
			}
		}
	}
	if (LBetaResidual != 0)
	{
		Betastart = LBeta - LBetaResidual;
		for (int numP = 0; numP < pTime; numP++)
		{
			pstart = numP * threadX;
			//mexPrintf("%d %d\n", pstart, Betastart);
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
				Filter << <block_cubic_residual, thread_cubic >> > (dev_Rcov, dev_R, dev_G, PInt, LP, pstart, Betastart, gstart, gend);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at pstart: %d Betastart: %d numg: %d\n",
						pstart, Betastart, numg);
					goto Error;
				}
			}
		}
	}

Error1:
	cudaFree(dev_R);
	cudaFree(dev_G);

	mexPrintf("backprojection\n");

	cudaStatus = cudaMalloc((void**)&dev_Size, 2 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Size cudaMalloc failed!\n");
		mexPrintf("dev_Size cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Size, Size, 2 * sizeof(double), cudaMemcpyHostToDevice);
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
	for (int BetaIndex = 0; BetaIndex < LBeta; BetaIndex++)
	{
		for (int numT = 0; numT < T_Time; numT++)
		{
			for (int numS = 0; numS < S_Time; numS++)
			{
				T_start = numT * threadX;
				S_start = numS * blockX;
				BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length,
					BetaScanRange[BetaIndex], dev_Pdomain, PInt, BetaScanInt, minP, maxP, BetaIndex, LP, LBeta, T_start, S_start, Rscan);
				
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at BetaIndex: %d\n", BetaIndex);
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
				BackProjection << <block_cubic_Bp, thread_cubic_Bp_residual >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length,
					BetaScanRange[BetaIndex], dev_Pdomain, PInt, BetaScanInt, minP, maxP, BetaIndex, LP, LBeta, T_start, S_start, Rscan);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at BetaIndex: %d\n", BetaIndex);
					goto Error;
				}
			}
			if (SlengthResidual != 0)
			{
				S_start = s_length - SlengthResidual;
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp_residual >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length,
					BetaScanRange[BetaIndex], dev_Pdomain, PInt, BetaScanInt, minP, maxP, BetaIndex, LP, LBeta, T_start, S_start, Rscan);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at BetaIndex: %d\n", BetaIndex);
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
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length,
					BetaScanRange[BetaIndex], dev_Pdomain, PInt, BetaScanInt, minP, maxP, BetaIndex, LP, LBeta, T_start, S_start, Rscan);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at BetaIndex: %d\n", BetaIndex);
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
	//for debug{
	//cudaStatus = cudaMemcpy(Rcov, dev_Rcov, LR * sizeof(float), cudaMemcpyDeviceToHost);
	//if (cudaStatus != cudaSuccess) {
	//	fprintf(stderr, "cudaMemcpy failed!\n");
	//	mexPrintf("cudaMemcpy dev_Display failed! %s\n", cudaGetErrorString(cudaStatus));
	//	goto Error;
	//}
	//}end debug

Error:
	cudaFree(dev_BetaScanRange);
	cudaFree(dev_Pdomain);
	cudaFree(dev_Rcov);
	cudaFree(dev_Display);
	cudaFree(dev_Size);

	mexPrintf("Exit FBP\n");
	return cudaStatus;
}