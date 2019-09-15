#include "BackprojectionFan2D.h"
// 2019/09/14 apply GPU acceleration
#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;
__device__ const double EPS = 1e-15;

__global__ void BackProjection(const float *dev_R, float *dev_Display, const double *dev_Size, const int T_length, const int S_length,
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
	double backweight = 0;

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
		Display_pBeta = ( P2 * dev_R[BetaIndex * LP + P1index] + P1 * dev_R[BetaIndex * LP + P2index] ) / PInt ;
	}

	double source_t = Center_t - Rscan * sin(Beta);      // define the source in matlab coordinate
	double source_s = Center_s + Rscan * cos(Beta);

	double Theta = atan(P_domain / Rscan);        // radian angle in s'-t coordinate plane 

	double Smax = 2 * Rscan;

	// define end detect point in matlab coordinate, Note that : 0 is the start
	double DetectPoint_tend = Center_t + Smax * sin(Theta) * cos(Beta) - (Rscan - Smax * cos(Theta)) * sin(Beta);
	double DetectPoint_send = Center_s + Smax * sin(Theta) * sin(Beta) + (Rscan - Smax * cos(Theta)) * cos(Beta);

	double T2S, S2T;
	if ((DetectPoint_tend - source_t) == 0)
		T2S = (DetectPoint_send - source_s) / (DetectPoint_tend - source_t + EPS);
	else
		T2S = (DetectPoint_send - source_s) / (DetectPoint_tend - source_t);
	if ((DetectPoint_send - source_s) == 0)
		S2T = (DetectPoint_tend - source_t) / (DetectPoint_send - source_s + EPS);
	else
		S2T = (DetectPoint_tend - source_t) / (DetectPoint_send - source_s);

	// limit the range of slope
	T2S = Maxlim(T2S); T2S = Minlim(T2S);
	S2T = Maxlim(S2T); S2T = Minlim(S2T);

	// actual Size
	double tlow = Tindex * Resolution_t, thigh = (Tindex+1) * Resolution_t, slow = Sindex * Resolution_s, shigh = (Sindex+1) * Resolution_s;

	//compute the first and last point in the ROI
	// using DetectPoint_end set up projection equation
	double tlow_s = source_s + (tlow - source_t) * T2S;
	double thigh_s = source_s + (thigh - source_t) * T2S;

	double slow_t = source_t + (slow - source_s) * S2T;
	double shigh_t = source_t + (shigh - source_s) * S2T;

	double T1 = 0, S1 = 0;
	double LengthinPixel =0;

	if (tlow_s >= 0 && tlow_s <= shigh)
	{
		T1 = tlow; S1 = tlow_s;
	}
	else if (slow_t >= 0 && slow_t <= thigh)
	{
		T1 = slow_t; S1 = slow;
	}
	else
	{
		//dev_Projection[threadid] = threadid;
		return;
	}

	LengthinPixel = 2 * Distance(T1, S1, DetectPoint_tend, DetectPoint_send);
	//if (LengthinROI == 0)
	//	return;
	backweight = LengthinPixel /*/ LengthinROI*/;   // no need for normalization
	thread_id = Sindex * T_length + Tindex;
	dev_Display[thread_id] += Display_pBeta * backweight;

}

// HeLTer function for using CUDA to add vectors in parallel.
cudaError_t FBPfan(float *Display, const float *R, const float *Pdomain, const float *BetaScanRange, const int LBeta,
	const int LP, const double *Size, const int t_length, const int s_length, const double Rscan)
{
	mexPrintf("Hello FBPpara!\n");
	float *dev_R = 0;
	float *dev_Pdomain = 0;
	double *dev_Size = 0;

	float *dev_Display = 0;
	float PInt = fabs(Pdomain[1] - Pdomain[0]);
	float BetaScanInt = fabs(BetaScanRange[1] - BetaScanRange[0]);

	float maxP = MAX(Pdomain[0], Pdomain[LP - 1]);
	float minP = MIN(Pdomain[0], Pdomain[LP - 1]);
	//mexPrintf("%lf %lf %lf %lf \n", maxGama, minGama, maxXigama, minXigama);

	const long LDisplay = t_length * s_length;
	const long LR = LP * LBeta;

	short thread_cubic_Bp_x = MIN(threadX, t_length);
	short block_cubic_Bp_x = MIN(blockX, s_length);

	const dim3 thread_cubic_Bp(thread_cubic_Bp_x, 1, 1);
	const dim3 block_cubic_Bp(block_cubic_Bp_x, 1, 1);

	dim3 thread_cubic_Bp_residual(1, 1, 1);  // initial
	dim3 block_cubic_Bp_residual(1, 1, 1);  // initial

	short TlengthResidual = t_length % threadX;
	short SlengthResidual = s_length % blockX;
	short T_Time = t_length / threadX;
	short S_Time = s_length / blockX;
	short T_start = 0;
	short S_start = 0;

	//mexPrintf("pTime: %d BetaTime: %d\n", pTime, BetaTime);

	if (TlengthResidual != 0)
	{
		thread_cubic_Bp_residual.x = TlengthResidual;
	}
	if (SlengthResidual != 0)
	{
		block_cubic_Bp_residual.x = SlengthResidual;
	}

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
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Pdomain, LP * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pdomain cudaMalloc failed!\n");
		mexPrintf("dev_Pdomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	//mexPrintf("copy data in CPU to GPU\n");

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_R, R, LR * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy R failed!\n");
		mexPrintf("cudaMemcpy R failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Pdomain, Pdomain, LP * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Pdomain failed!\n");
		mexPrintf("cudaMemcpy Pdomain failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	mexPrintf("start parallel computation\n");

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
				BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_R, dev_Display, dev_Size, t_length, s_length,
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
				BackProjection << <block_cubic_Bp, thread_cubic_Bp_residual >> > (dev_R, dev_Display, dev_Size, t_length, s_length,
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
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp_residual >> > (dev_R, dev_Display, dev_Size, t_length, s_length,
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
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp >> > (dev_R, dev_Display, dev_Size, t_length, s_length,
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
	cudaFree(dev_Pdomain);
	cudaFree(dev_R);
	cudaFree(dev_Display);
	cudaFree(dev_Size);

	mexPrintf("Exit FBP\n");
	return cudaStatus;
}