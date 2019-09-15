#include "IterativeFDKerrCompute.h"
#include <cstdlib>
// 2018/08/07
cudaError_t IterativeFDKerrCompute(const float *Pic, float *ErrorSlicecuda, const float *BetaScanRange, const float *Udomain,
	const double *Pdomain, const double *G, const double *Size, const int t_length, const int s_length, const int z_length,
	const double Center_t, const double Center_s, const double Center_z, const int LBeta, const int LP, const int LU,
	const int LG, const double Distance, const double *resolution2, const int t, const float BetaScanInt, const float dU,
	const double PInt);

__device__ const double PI = 3.141592653589793;

__global__ void IterativeFDKerr(const float *dev_Pic, float *dev_ErrorSlicecuda, const float U, const double *dev_Pdomain,
	 const double *dev_G, const double *dev_Size, const double *dev_resolution2, const double Center_t,
	const double Center_s, const double Center_z, const int t_length, const int s_length, 
	const int z_length, const float Beta, const int numB, const int sstart, const int zstart, const double Distance, 
	const int LBeta, const int LP, const int LU, const int LG, const int t, const float BetaScanInt, const float dU, const double PInt)
	
{
	const unsigned int s = threadIdx.x + sstart;
	const unsigned int z = blockIdx.x + zstart;
	 
	const unsigned long threadid = numB * s_length * z_length + z * s_length + s;
	//const unsigned long threadid = z * s_length * LBeta + s * LBeta + numB;

	float image_t = 0, image_s = 0, image_z = 0, dect_t = 0, dect_s = 0, dect_z = 0, LengthRatio = 0, T_deriv = 0, S_deriv = 0,
		Z_deriv = 0, X_deriv = 0, Y_deriv = 0;
	double P = 0 ;

	int filter_index = 0, X_deriv_index = 0, Y_deriv_index = 0, Z_deriv_index = 0, Pic_index1 = 0, Pic_index2 = 0;

	// image pixel in ground coordinate
	image_t = (t + 0.5) * dev_resolution2[0] - Center_t;
	image_s = (s + 0.5) * dev_resolution2[1] - Center_s;
	image_z = (z + 0.5) * dev_resolution2[2] - Center_z;

	dect_t = image_t * cos(Beta) + image_s * sin(Beta);          // in rotate coordinate
	dect_s = -image_t * sin(Beta) + image_s * cos(Beta);
	dect_z = image_z;

	LengthRatio = Distance / (Distance - dect_s);
	for (int p = 0; p<LP; p++)
	{
		P = dev_Pdomain[p];
		filter_index = floor((dect_t * LengthRatio - P) / PInt) + LP;
		if (filter_index >= 0 && filter_index < LG)
		{
		/*	for (int u = 0; u<LU; u++)
			{*/
				T_deriv = P * (Distance - U) / Distance;                         // in rotate coordinate
				S_deriv = U;
				Z_deriv = dect_z + (dect_s - U) * dect_z / (Distance - dect_s);

				X_deriv = T_deriv * cos(Beta) - S_deriv * sin(Beta);          // in ground coordinate
				Y_deriv = T_deriv * sin(Beta) + S_deriv * cos(Beta);

				X_deriv_index = floor((X_deriv + Center_t) / dev_resolution2[0]);
				Y_deriv_index = floor((Y_deriv + Center_s) / dev_resolution2[1]);
				Z_deriv_index = floor((Z_deriv + Center_z) / dev_resolution2[2]);
				
				
				if (X_deriv_index >= 0 && X_deriv_index < t_length  &&  Y_deriv_index >= 0
					&& Y_deriv_index < s_length && Z_deriv_index >= 0 && Z_deriv_index < z_length)
				{
					Pic_index1 = Z_deriv_index * (t_length * s_length) + Y_deriv_index * t_length + X_deriv_index;
					Pic_index2 = z * (t_length * s_length) + Y_deriv_index * t_length + X_deriv_index;

					//dev_ErrorSlicecuda[threadid] = BetaScanInt ;
					dev_ErrorSlicecuda[threadid] += 0.5 * pow2(LengthRatio) * (dev_Pic[Pic_index1]
						- dev_Pic[Pic_index2]) * sqrt(1 + pow2(P / Distance))
						* dev_G[filter_index] * Distance / sqrt(pow2(Distance) + pow2(P)) * dU * PInt * BetaScanInt;
				}//if-xyzindex
			//}// for-u
		}//if-filter
	}//for-p
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t IterativeFDKerrCompute(const float *Pic, float *ErrorSlicecuda, const float *BetaScanRange, const float *Udomain,
	const double *Pdomain, const double *G, const double *Size, const int t_length, const int s_length, const int z_length,
	const double Center_t, const double Center_s, const double Center_z, const int LBeta, const int LP, const int LU,
	const int LG, const double Distance, const double *resolution2, const int t, const float BetaScanInt, const float dU,
	const double PInt)
{
	mexPrintf("Start CPU-GPU interface!\n");
	
	float *dev_Pic = 0, *dev_Udomain = 0;
	float *dev_ErrorSlicecuda = 0;
	double *dev_Pdomain = 0;
	double *dev_resolution2 = 0, *dev_Size = 0, *dev_G = 0 ;

	int threadcubic_x = MIN(threadX, s_length);
	int blockcubic_x = MIN(blockX, z_length);
	int slresidual = s_length % threadX;
	int zlresidual = z_length % blockX;
	int sTime = s_length / threadX;
	int zTime = z_length / blockX;
	int sstart = 0, zstart = 0;
	float beta = 0, U = 0; double P = 0;
	const dim3 thread_cubic(threadcubic_x,1,1);     // num of threads in each block depends on s_length
	const dim3 block_cubic(blockcubic_x, 1, 1);      // num of blocks in each grid depends on s_length
	dim3 thread_cubic_residual(1, 1, 1);    // initial, size of dim3 depends on slresidual/zlresidual
	dim3 block_cubic_residual(1, 1, 1);  
	
	float *Zeros = new float[LBeta * s_length * z_length];      // for initialize the output
	memset(Zeros, 0, sizeof(float));

	if (slresidual != 0)
	{
		thread_cubic_residual.x = slresidual;
	}
	if (zlresidual != 0)
	{
		block_cubic_residual.x = zlresidual;
	}

    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		mexPrintf("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

	mexPrintf("Start to call for GPU space.\n");

    // Allocate GPU buffers for three vectors (6 input, 1 output)    .
    cudaStatus = cudaMalloc((void**)&dev_Pic, t_length * s_length * z_length * sizeof(float));
    if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pic cudaMalloc failed!");
		mexPrintf("dev_Pic cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_Udomain, LU * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_Udomain cudaMalloc failed!");
		mexPrintf("dev_Udomain cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_Pdomain, LP * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_Pdomain cudaMalloc failed!");
		mexPrintf("dev_Pdomain cudaMalloc failed!");
        goto Error;
    }

	cudaStatus = cudaMalloc((void**)&dev_G, LG * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_G cudaMalloc failed!");
		mexPrintf("dev_G cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_resolution2, 3 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution2 cudaMalloc failed!");
		mexPrintf("dev_resolution2 cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Size, 3 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Size cudaMalloc failed!");
		mexPrintf( "dev_Size cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_ErrorSlicecuda, LBeta * s_length * z_length * sizeof(float));    // output
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_ErrorSlice cuda cudaMalloc failed!");
		mexPrintf("dev_ErrorSlice cuda cudaMalloc failed!");
		goto Error;
	}

	mexPrintf("Start to copy data from CPU to GPU space.\n");

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_Pic, Pic, t_length * s_length * z_length * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_Pic cudaMemcpy failed!");
		mexPrintf("dev_Pic cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_Udomain, Udomain, LU * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_Udomain cudaMemcpy failed!");
		mexPrintf("dev_Udomain cudaMemcpy failed!");
        goto Error;
    }

	cudaStatus = cudaMemcpy(dev_Pdomain, Pdomain, LP * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pdomain cudaMemcpy failed!");
		mexPrintf("dev_Pdomain cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_G, G, LG * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_G cudaMemcpy failed!");
		mexPrintf("dev_G cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_resolution2, resolution2, 3 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution2 cudaMemcpy failed!");
		mexPrintf("dev_resolution2 cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Size, Size, 3 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Size cudaMemcpy failed!");
		mexPrintf("dev_Size cudaMemcpy failed!");
		goto Error;
	}

	// initialize the output
	cudaStatus = cudaMemcpy(dev_ErrorSlicecuda, Zeros, LBeta * s_length * z_length * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_ErrorSlicecuda initialization failed!");
		mexPrintf("dev_ErrorSlicecuda initialization failed!");
		goto Error;
	}

	mexPrintf("start to compute discrepancy occur in each point.\n");

    // Launch a kernel on the GPU with one thread for each element.
	for (int numB = 0; numB < LBeta; numB++)            // rotation angle
	{
		beta = BetaScanRange[numB];
		//mexPrintf("numbeta: %d\n", numB);
		for (int numU = 0; numU < LU; numU++)      // row on the detector
		{
			U = Udomain[numU];
			//for (int numP = 0; numP < LP; numP++)      // row on the detector
			//{
			//	P = Pdomain[numP];
		
				for (int nums = 0; nums < sTime; nums++)
				{
					for (int numz = 0; numz < zTime; numz++)
					{
						sstart = nums * threadX;
						zstart = numz * blockX;
						//mexPrintf("%d %d\n", sstart, zstart);

						IterativeFDKerr << <block_cubic, thread_cubic >> >(dev_Pic, dev_ErrorSlicecuda, U, dev_Pdomain, dev_G,
							dev_Size, dev_resolution2, Center_t, Center_s, Center_z, t_length, s_length, z_length, beta, numB, sstart,
							zstart, Distance, LBeta, LP, LU, LG, t, BetaScanInt, dU, PInt);
					}
				}

				if (slresidual != 0)
				{
					sstart = s_length - slresidual;
					if (zlresidual != 0)
					{
						zstart = z_length - zlresidual;
						//("%d %d\n", sstart, zstart);
						IterativeFDKerr << <block_cubic, thread_cubic >> >(dev_Pic, dev_ErrorSlicecuda, U, dev_Pdomain, dev_G,
							dev_Size, dev_resolution2, Center_t, Center_s, Center_z, t_length, s_length, z_length, beta, numB, sstart,
							zstart, Distance, LBeta, LP, LU, LG, t, BetaScanInt, dU, PInt);
					}

					for (int numz = 0; numz < zTime; numz++)
					{
						zstart = numz * blockX;
						//("%d %d\n", sstart, zstart);
						IterativeFDKerr << <block_cubic, thread_cubic >> >(dev_Pic, dev_ErrorSlicecuda, U, dev_Pdomain, dev_G,
							dev_Size, dev_resolution2, Center_t, Center_s, Center_z, t_length, s_length, z_length, beta, numB, sstart,
							zstart, Distance, LBeta, LP, LU, LG, t, BetaScanInt, dU, PInt);
					}
				}
				if (zlresidual != 0)
				{
					zstart = z_length - zlresidual;
					for (int nums = 0; nums < sTime; nums++)
					{
						sstart = nums * threadX;
						//mexPrintf("%d %d\n", sstart, zstart);
						IterativeFDKerr << <block_cubic, thread_cubic >> >(dev_Pic, dev_ErrorSlicecuda, U, dev_Pdomain, dev_G,
							dev_Size, dev_resolution2, Center_t, Center_s, Center_z, t_length, s_length, z_length, beta, numB, sstart,
							zstart, Distance, LBeta, LP, LU, LG, t, BetaScanInt, dU, PInt);
					}
				}
				
			//}//for-nump
		}//for-numU		
	}//for-numB

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		mexPrintf("addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		mexPrintf("cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(ErrorSlicecuda, dev_ErrorSlicecuda, LBeta * s_length * z_length * sizeof(float), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
		mexPrintf("ErrorSlicecuda cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_ErrorSlicecuda);
    cudaFree(dev_Pic);
    cudaFree(dev_Udomain);
	cudaFree(dev_Pdomain);
	cudaFree(dev_G);
	cudaFree(dev_Size);
	cudaFree(dev_resolution2);

    return cudaStatus;
}
