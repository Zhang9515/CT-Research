#include "FBPparallel.h"
// 2019/03/28 apply GPU acceleration
#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;

__global__ void GFunction(float *dev_G, const float TInt, const int tstart, const int LT)
{

	const unsigned int Tindex = tstart + threadIdx.x;

	// S_L filter 
	if (Tindex == 0)
		dev_G[LT - 1] = 1.0 / pow2(PI*TInt);
	else
	{
		dev_G[(LT - 1) - Tindex] = -1.0 / pow2(TInt * PI) / (4 * pow2(Tindex) - 1);
		dev_G[(LT - 1) + Tindex] = -1.0 / pow2(TInt * PI) / (4 * pow2(Tindex) - 1);
	}
}

__global__ void Filter(float *dev_Rcov, float *dev_R, const float *dev_G, const float TInt, 
	const int LT, const int tstart, const int Thetastart, const int gstart, const int gend)
{
	const unsigned int Tindex = threadIdx.x + tstart;
	const unsigned int thetaIndex = blockIdx.x + Thetastart;
	const unsigned long base_id = thetaIndex * LT;
	const unsigned long thread_id = base_id + Tindex;

	double Rcovsum = 0;
	// convolution

	for (int g = gstart; g < gend; g++)
	{
		//convolution
		Rcovsum += dev_R[base_id + g] * dev_G[(LT - 1) + Tindex - g];
	}

	dev_Rcov[thread_id] += TInt * Rcovsum;
}

__global__ void BackProjection(const float *dev_Rcov, float *dev_Display, const double *dev_Size, const int T_length, const int S_length,
	const float Theta, const float *dev_Tdomain, const float TInt, const float ThetaScanInt, const float mint, 
	const float maxt, const int thetaIndex, const int LT, const int LTheta, const int T_start, const int S_start)
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

	double image_t, image_s, dect_t;
	unsigned int tN1index, tN2index;
	double t_domain1, t_domain2;
	double t1, t2;
	double Display_pTheta = 0;
	double t;
	unsigned long thread_id;

	// this is a little different from code on MATLAB
	image_t = (Tindex + 0.5) * Resolution_t - Center_t;    // center of the image pixel in ground coordinate
	image_s = (Sindex + 0.5) * Resolution_s - Center_s;

	t = image_t * cos(Theta) + image_s * sin(Theta);   // rotate in ground coordinate
	//dect_t = -image_t * sin(Theta) + image_s * cos(Theta);   // rotate in ground coordinate
	
	if ((t >= mint) && (t < maxt))
	{
		tN1index = floor(fabs(t - dev_Tdomain[0]) / TInt);
		tN2index = tN1index + 1;

		t_domain1 = dev_Tdomain[tN1index]; t_domain2 = dev_Tdomain[tN2index];
	
		//bilinear interpolation
		t1 = fabs(t - t_domain1); t2 = fabs(t_domain2 - t);

		//Display_pTheta = 1;
		Display_pTheta = ( t2 * dev_Rcov[thetaIndex * LT + tN1index] + t1 * dev_Rcov[thetaIndex * LT + tN2index] )
			/ TInt * ThetaScanInt;
	}
	thread_id = Sindex * T_length + Tindex;
	dev_Display[thread_id] += Display_pTheta;

}

// HeLTer function for using CUDA to add vectors in parallel.
cudaError_t FBPpara(float *Display, const float *R, const float *Tdomain, const float *ThetaScanRange, const int LTheta,
	const int LT, const double *Size, const int t_length, const int s_length)
{
	mexPrintf("Hello FBPpara!\n");
	float *dev_R = 0, *dev_Rcov = 0;
	float *dev_G = 0;
	float *dev_ThetaScanRange = 0, *dev_Tdomain = 0;
	double *dev_Size = 0;

	float *dev_Display = 0;
	float TInt = fabs(Tdomain[1] - Tdomain[0]);
	float ThetaScanInt = fabs(ThetaScanRange[1] - ThetaScanRange[0]);

	float maxt = MAX(Tdomain[0], Tdomain[LT - 1]);
	float mint = MIN(Tdomain[0], Tdomain[LT - 1]);
	//mexPrintf("%lf %lf %lf %lf \n", maxGama, minGama, maxXigama, minXigama);

	const long LDisplay = t_length * s_length;
	const long LR = LT * LTheta;
	const int LFilter = 2 * LT - 1;

	short thread_cubic_x = MIN(threadX, LT);
	short block_cubic_x = MIN(blockX, LTheta);
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

	short LTResidual = LT % threadX;
	short LThetaResidual = LTheta % blockX;
	short tTime = LT / threadX;
	short ThetaTime = LTheta / blockX;
	short tstart = 0;
	short Thetastart = 0;
	short TlengthResidual = t_length % threadX;
	short SlengthResidual = s_length % blockX;
	short T_Time = t_length / threadX;
	short S_Time = s_length / blockX;
	short T_start = 0;
	short S_start = 0;
	short gstart = 0, gend = 0;
	short gtime = 1;
	//mexPrintf("tTime: %d ThetaTime: %d\n", tTime, ThetaTime);

	if (LTResidual != 0)
	{
		thread_cubic_residual.x = LTResidual;
	}
	if (LThetaResidual != 0)
	{
		block_cubic_residual.x = LThetaResidual;
	}
	if (TlengthResidual != 0)
	{
		thread_cubic_Bp_residual.x = TlengthResidual;
	}
	if (SlengthResidual != 0)
	{
		block_cubic_Bp_residual.x = SlengthResidual;
	}

	mexPrintf("thread_cubic_x: %d block_cubic_x: %d LTResidual: %d LThetaResidual: %d\n",
		thread_cubic_x, block_cubic_x, LTResidual, LThetaResidual);
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

	cudaStatus = cudaMalloc((void**)&dev_Tdomain, LT * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Tdomain cudaMalloc failed!\n");
		mexPrintf("dev_Tdomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_ThetaScanRange, LTheta * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_ThetaScanRange cudaMalloc failed!\n");
		mexPrintf("dev_ThetaScanRange cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
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

	cudaStatus = cudaMemcpy(dev_Tdomain, Tdomain, LT * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Tdomain failed!\n");
		mexPrintf("cudaMemcpy Tdomain failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_ThetaScanRange, ThetaScanRange, LTheta * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy ThetaScanRange failed!\n");
		mexPrintf("cudaMemcpy ThetaScanRange failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	mexPrintf("start parallel computation\n");

	mexPrintf("Preweight and Filtering\n");
	// Preweight and filter 
	// to be adapted to GPU, we limited the number of threads in each grid (threadX * blockX * LTheta)
	// generate filter
	mexPrintf("Generate filter\n");
	cudaStatus = cudaMalloc((void**)&dev_G, LFilter * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_G cudaMalloc failed!\n");
		mexPrintf("dev_G cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error1;
	}
	cudaMemset(dev_G, 0, sizeof(float));

	for (int numT = 0; numT < tTime; numT++)
	{
		tstart = numT * threadX;
		GFunction << <1, thread_cubic_x >> > (dev_G, TInt, tstart, LT);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "generate filter failed: %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("generate filter launch failed %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("Error happens at numT: %d\n", numT);
			goto Error;
		}
	}
	if (LTResidual != 0)
	{
		tstart = LT - LTResidual;
		GFunction << <1, thread_cubic_residual.x >> > (dev_G, TInt, tstart, LT);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "generate filter failed: %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("generate filter launch failed %s\n", cudaGetErrorString(cudaStatus));
			mexPrintf("Error happens at LTResidual: %d\n", LTResidual);
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
	gtime = ceil(1.0 * LT / Filterlengthlimit);
	mexPrintf("gtime: %d\n", gtime);

	for (int numT = 0; numT < tTime; numT++)
	{
		for (int numTheta = 0; numTheta < ThetaTime; numTheta++)
		{
			tstart = numT * threadX;
			Thetastart = numTheta * blockX;
			//mexPrintf("%d %d\n", tstart, Thetastart);
			for (int numg = 0; numg < gtime; numg++)
			{
				gstart = int(numg * Filterlengthlimit);
				if (numg != (gtime - 1))
				{
					gend = int((numg + 1)* Filterlengthlimit);
				}
				else
				{
					gend = LT;
				}
				/*mexPrintf("gstart: %d\n", gstart);
				mexPrintf("gend: %d\n", gend);*/
				Filter << <block_cubic, thread_cubic >> > (dev_Rcov, dev_R, dev_G, TInt, LT, tstart, Thetastart, gstart, gend);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at tstart: %d Thetastart: %d numg: %d\n",
						tstart, Thetastart, numg);
					goto Error;
				}
			}
		}
	}

	if (LTResidual != 0)
	{
		tstart = LT - LTResidual;
		if (LThetaResidual != 0)
		{
			Thetastart = LTheta - LThetaResidual;
			//("%d %d\n", tstart, Thetastart);
			for (int numg = 0; numg < gtime; numg++)
			{
				gstart = int(numg * Filterlengthlimit);
				if (numg != (gtime - 1))
				{
					gend = int((numg + 1)* Filterlengthlimit);
				}
				else
				{
					gend = LT;
				}
				/*mexPrintf("gstart: %d\n", gstart);
				mexPrintf("gend: %d\n", gend);*/
				Filter << <block_cubic_residual, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_G, TInt, LT, tstart, Thetastart, gstart, gend);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at tstart: %d Thetastart: %d numg: %d\n",
						tstart, Thetastart, numg);
					goto Error;
				}
			}
		}

		for (int numTheta = 0; numTheta < ThetaTime; numTheta++)
		{
			Thetastart = numTheta * blockX;
			//("%d %d\n", tstart, Thetastart);
			for (int numg = 0; numg < gtime; numg++)
			{
				gstart = int(numg * Filterlengthlimit);
				if (numg != (gtime - 1))
				{
					gend = int((numg + 1)* Filterlengthlimit);
				}
				else
				{
					gend = LT;
				}
				/*mexPrintf("gstart: %d\n", gstart);
				mexPrintf("gend: %d\n", gend);*/
				Filter << <block_cubic, thread_cubic_residual >> > (dev_Rcov, dev_R, dev_G, TInt, LT, tstart, Thetastart, gstart, gend);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at tstart: %d Thetastart: %d numg: %d\n",
						tstart, Thetastart, numg);
					goto Error;
				}
			}
		}
	}
	if (LThetaResidual != 0)
	{
		Thetastart = LTheta - LThetaResidual;
		for (int numT = 0; numT < tTime; numT++)
		{
			tstart = numT * threadX;
			//mexPrintf("%d %d\n", tstart, Thetastart);
			for (int numg = 0; numg < gtime; numg++)
			{
				gstart = int(numg * Filterlengthlimit);
				if (numg != (gtime - 1))
				{
					gend = int((numg + 1)* Filterlengthlimit);
				}
				else
				{
					gend = LT;
				}
				/*mexPrintf("gstart: %d\n", gstart);
				mexPrintf("gend: %d\n", gend);*/
				Filter << <block_cubic_residual, thread_cubic >> > (dev_Rcov, dev_R, dev_G, TInt, LT, tstart, Thetastart, gstart, gend);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Filter launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Filter launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at tstart: %d Thetastart: %d numg: %d\n",
						tstart, Thetastart, numg);
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
	for (int thetaIndex = 0; thetaIndex < LTheta; thetaIndex++)
	{
		for (int numT = 0; numT < T_Time; numT++)
		{
			for (int numS = 0; numS < S_Time; numS++)
			{
				T_start = numT * threadX;
				S_start = numS * blockX;
				BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length,
					ThetaScanRange[thetaIndex], dev_Tdomain, TInt, ThetaScanInt, mint, maxt, thetaIndex, LT, LTheta, T_start, S_start);
				
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at thetaIndex: %d\n", thetaIndex);
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
					ThetaScanRange[thetaIndex], dev_Tdomain, TInt, ThetaScanInt, mint, maxt, thetaIndex, LT, LTheta, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at thetaIndex: %d\n", thetaIndex);
					goto Error;
				}
			}
			if (SlengthResidual != 0)
			{
				S_start = s_length - SlengthResidual;
				BackProjection << <block_cubic_Bp_residual, thread_cubic_Bp_residual >> > (dev_Rcov, dev_Display, dev_Size, t_length, s_length,
					ThetaScanRange[thetaIndex], dev_Tdomain, TInt, ThetaScanInt, mint, maxt, thetaIndex, LT, LTheta, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at thetaIndex: %d\n", thetaIndex);
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
					ThetaScanRange[thetaIndex], dev_Tdomain, TInt, ThetaScanInt, mint, maxt, thetaIndex, LT, LTheta, T_start, S_start);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("BackProjection launch failed %s\n", cudaGetErrorString(cudaStatus));
					mexPrintf("Error happens at thetaIndex: %d\n", thetaIndex);
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
	cudaFree(dev_ThetaScanRange);
	cudaFree(dev_Tdomain);
	cudaFree(dev_Rcov);
	cudaFree(dev_Display);
	cudaFree(dev_Size);

	mexPrintf("Exit FBP\n");
	return cudaStatus;
}