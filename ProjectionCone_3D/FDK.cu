#include "FDK.h"
// 2018/04/20 apply GPU acceleration
#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;

__global__ void GFunction(float *dev_G, const float PInt, const int LP)
{

	const unsigned int Pindex = threadIdx.x;

	// S_L filter 
	if (Pindex == 0)
		dev_G[LP - 1] = 1.0 / pow2(PI*PInt);
	else
	{
		dev_G[(LP - 1) - Pindex] = -1.0 / pow2(PInt * PI) / (4 * pow2(Pindex) - 1);
		dev_G[(LP - 1) + Pindex] = -1.0 / pow2(PInt * PI) / (4 * pow2(Pindex) - 1);
	}
}

__global__ void PreWeightFiltration(float *dev_Rcov, float *dev_R, const float *dev_G, const float *dev_Pdomain,
	const float *dev_Xigamadomain, const double Distance, const float PInt, const int LP, const int LXigama, const int Pstart,
	const int Xigamastart, const int Betaindex)
{

	//const unsigned int Gamaindex = blockIdx.x * blockDim.x + threadIdx.x;
	//const unsigned int Xigamaindex = blockIdx.y * blockDim.y + threadIdx.y;
	//const unsigned int Betaindex = blockIdx.z * blockDim.z + threadIdx.z;
	//const unsigned long thread_id = Betaindex * ( gridDim.x * gridDim.y * blockDim.x * blockDim.y ) 
	//	+ Xigamaindex * ( gridDim.x * blockDim.x ) + Gamaindex ;

	const unsigned int Pindex = threadIdx.x + Pstart;
	const unsigned int Xigamaindex = blockIdx.x + Xigamastart;
	const unsigned long thread_id = Betaindex * (LXigama * LP) + Xigamaindex * LP + Pindex;

	float P = dev_Pdomain[Pindex];
	float Xigama = dev_Xigamadomain[Xigamaindex];

	double Proportion = Distance / sqrt(pow2(Distance) + pow2(P) + pow2(Xigama));

	dev_R[thread_id] = dev_R[thread_id] * Proportion;    // directly cover the input

	__syncthreads();
	double Rcovsum = 0;
	// convolution
	for (int g = 0; g < LP; g++)
	{
		//convolution
		Rcovsum += dev_R[Betaindex * (LXigama * LP) + Xigamaindex * LP + g] * dev_G[(LP - 1) + Pindex - g];
	}
	//__syncthreads();
	dev_Rcov[thread_id] = PInt * Rcovsum;
}

__global__ void BackProjection(const float *dev_Rcov, float *dev_Display, const double *dev_Size,
	const int t_length, const int s_length, const int z_length, const float Beta, const double Distance, 
	const float *dev_Pdomain, const float *dev_Xigamadomain, const float PInt, const float XigamaInt, 
	const float BetaScanInt, const float minP, const float maxP, const float minXigama, const float maxXigama, 
	const int betaIndex, const int LP, const int LXigama)
{
	const unsigned int Tindex = threadIdx.x;
	const unsigned int Sindex = blockIdx.x;
	const unsigned int Zindex = blockIdx.y;
	const unsigned long thread_id = Zindex * (gridDim.x * blockDim.x)
		+ Sindex * blockDim.x + Tindex;
	// initialize

	const double Resolution_t = 1.0 * dev_Size[0] / t_length;
	const double Resolution_s = 1.0 * dev_Size[1] / s_length;
	const double Resolution_z = 1.0 * dev_Size[2] / z_length;

	// rotation center
	double Center_t = dev_Size[0] / 2;
	double Center_s = dev_Size[1] / 2;
	double Center_z = dev_Size[2] / 2;

	// this is a little different from code on MATLAB
	double image_t = (Tindex + 0.5) * Resolution_t - Center_t;  float image_s = (Sindex + 0.5) * Resolution_s - Center_s; float image_z = (Zindex + 0.5) * Resolution_z - Center_z;           // image pixel in ground coordinate

	// rotate in ground coordinate
	double dect_t = image_t * cos(Beta) + image_s * sin(Beta);
	double dect_s = -image_t * sin(Beta) + image_s * cos(Beta);
	double dect_z = image_z;

	// define the projection position on the detector
	double LengthRatio = Distance / (Distance - dect_s);
	double Xigama = dect_z * LengthRatio;
	double P = dect_t * LengthRatio;

	unsigned int XigamaN1index = 0, XigamaN2index = 0, PN1index = 0, PN2index = 0;
	double P_domain1 = 0, P_domain2 = 0, Xigama_domain1 = 0, Xigama_domain2 = 0;
	double Xig1 = 0, Xig2 = 0, P1 = 0, P2 = 0;
	double Display_pBeta = 0;

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

		double Weight = sqrt(1 + 0 * pow2(dect_z) / (pow2(Distance + dect_s) + pow2(dect_t)));
		//Display_pBeta = 1;
		Display_pBeta = (Xig2 * P2 * dev_Rcov[betaIndex * LP * LXigama + XigamaN1index * LP + PN1index]
			+ Xig1 * P2 * dev_Rcov[betaIndex * LP * LXigama + XigamaN2index * LP + PN1index] + Xig2 * P1 * dev_Rcov[betaIndex * LP * LXigama + XigamaN1index * LP + PN2index]
			+ Xig1 * P1 * dev_Rcov[betaIndex * LP * LXigama + XigamaN2index * LP  + PN2index]) / (PInt * XigamaInt) * pow2(LengthRatio) * BetaScanInt * Weight;
	}

	dev_Display[thread_id] += Display_pBeta;
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

	int thread_cubic_x = MIN(threadX, LP);
	int block_cubic_x = MIN(blockX, LXigama);

	const dim3 thread_cubic(thread_cubic_x, 1, 1);
	const dim3 block_cubic(block_cubic_x, 1, 1);

	dim3 thread_cubic_residual(1, 1, 1);  // initial
	dim3 block_cubic_residual(1, 1, 1);  // initial

	int LPResidual = LP % threadX;
	int LXigamaResidual = LXigama % blockX;
	int PTime = LP / threadX;
	int XigamaTime = LXigama / blockX;
	int Pstart = 0;
	int Xigamastart = 0;
	float Beta = 0;

	if (LPResidual != 0)
	{
		thread_cubic_residual.x = LPResidual;
	}
	if (LXigamaResidual != 0)
	{
		block_cubic_residual.x = LXigamaResidual;
	}

	const dim3 thread_cubic_Bp(t_length, 1, 1);
	const dim3 block_cubic_Bp(s_length, z_length, 1);

	cudaError_t cudaStatus;

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

	cudaStatus = cudaMalloc((void**)&dev_G, LFilter * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_G cudaMalloc failed!\n");
		mexPrintf("dev_G cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error1;
	}

	cudaStatus = cudaMalloc((void**)&dev_Rcov, LR * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_BetaScanRange cudaMalloc failed!\n");
		mexPrintf("dev_BetaScanRange cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Xigamadomain, LXigama * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Xigamadomain cudaMalloc failed!\n");
		mexPrintf("dev_Xigamadomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
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

	cudaStatus = cudaMemcpy(dev_BetaScanRange, BetaScanRange, LBeta * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy BetaScanRange failed!\n");
		mexPrintf("cudaMemcpy v failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Xigamadomain, Xigamadomain, LXigama * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Xigamadomain failed!\n");
		mexPrintf("cudaMemcpy Xigamadomain failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	mexPrintf("start parallel computation\n");
	mexPrintf("Filter\n");
	// Generate Filter
	GFunction << <1, LP >> >(dev_G, PInt, LP);

	mexPrintf("Preweight and filtration\n");
	// Preweight and filtration
	// to be adapted to GPU, we limited the number of threads in each grid (threadX * blockX * LBeta)
	for (int numB = 0; numB < LBeta; numB++)
	{
		Beta = BetaScanRange[numB];
		for (int numP = 0; numP < PTime; numP++)
		{
			for (int numX = 0; numX < XigamaTime; numX++)
			{
				Pstart = numP * threadX;
				Xigamastart = numX * blockX;
				//mexPrintf("%d %d\n", Pstart, Xigamastart);
				PreWeightFiltration << <block_cubic, thread_cubic >> > (dev_Rcov, dev_R, dev_G, dev_Pdomain,
					dev_Xigamadomain, Distance, PInt, LP, LXigama, Pstart, Xigamastart, numB);
			}
		}

		if (LPResidual != 0)
		{
			Pstart = LP - LPResidual;
			if (LXigamaResidual != 0)
			{
				Xigamastart = LXigama - LXigamaResidual;
				//("%d %d\n", Pstart, Xigamastart);
				PreWeightFiltration << <block_cubic_residual, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_G, dev_Pdomain,
					dev_Xigamadomain, Distance, PInt, LP, LXigama, Pstart, Xigamastart, numB);
			}

			for (int numX = 0; numX < XigamaTime; numX++)
			{
				Xigamastart = numX * blockX;
				//("%d %d\n", Pstart, Xigamastart);
				PreWeightFiltration << <block_cubic, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_G, dev_Pdomain,
					dev_Xigamadomain, Distance, PInt, LP, LXigama, Pstart, Xigamastart, numB);
			}
		}
		if (LXigamaResidual != 0)
		{
			Xigamastart = LXigama - LXigamaResidual;
			for (int numP = 0; numP < PTime; numP++)
			{
				Pstart = numP * threadX;
				//mexPrintf("%d %d\n", Pstart, Xigamastart);
				PreWeightFiltration << <block_cubic_residual, thread_cubic >> > (dev_Rcov, dev_R, dev_G, dev_Pdomain,
					dev_Xigamadomain, Distance, PInt, LP, LXigama, Pstart, Xigamastart, numB);
			}
		}
	}
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "PreWeight and Filtration failed: %s\n", cudaGetErrorString(cudaStatus));
		mexPrintf("PreWeight and Filtration failed\n");
		goto Error1;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error1;
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
		BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length, z_length,
			BetaScanRange[betaIndex], Distance, dev_Pdomain, dev_Xigamadomain, PInt, XigamaInt, BetaScanInt, minP, maxP, 
			minXigama, maxXigama, betaIndex, LP, LXigama);
	}

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
		mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
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