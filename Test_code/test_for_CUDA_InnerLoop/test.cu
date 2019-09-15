#include "test.h"
// 2018/04/20 apply GPU acceleration
#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;

__global__ void BackProjection(const float *dev_R, float *dev_Display, bool * dev_signal)
{

	const unsigned int Tindex = threadIdx.x;
	const unsigned int Bindex = blockIdx.x;
	int index = Bindex * 256 + Tindex;
	__syncthreads();
	//for (int num = 0;num<256*16;num++)
	//{
	//	if (num == 0)
	//	{
			while (1) {
				if (!dev_signal[Tindex]) {
					dev_signal[Tindex] = true;
					dev_Display[Tindex] += /*dev_R[index]*/Bindex;
					dev_signal[Tindex] = false;
					__threadfence();
					break;
				}
			}
				

	//	}		
	//}

	//__threadfence();
	//dev_Display[Tindex] = 3;

}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t FDKpro(float *Display, const float *R)
{
	float* dev_Display = 0; bool *dev_signal = 0; float* dev_R = 0;
	int LR = 16 * 256; int LD = 256;

	const dim3 thread_cubic(256, 1, 1);
	const dim3 block_cubic(16, 1, 1);

	cudaError_t cudaStatus;

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

	cudaStatus = cudaMalloc((void**)&dev_Display, LD * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pdomain cudaMalloc failed!\n");
		mexPrintf("dev_Pdomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	cudaMemset(dev_Display, 0, sizeof(float));

	cudaStatus = cudaMalloc((void**)&dev_signal, LD * sizeof(bool));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pdomain cudaMalloc failed!\n");
		mexPrintf("dev_Pdomain cudaMalloc failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	cudaMemset(dev_signal, false, sizeof(bool));
	//mexPrintf("copy data in CPU to GPU\n");

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_R, R, LR * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy R failed!\n");
		mexPrintf("cudaMemcpy R failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	//Backprojection
	BackProjection << <block_cubic, thread_cubic >> > (dev_R, dev_Display, dev_signal);

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
	cudaStatus = cudaMemcpy(Display, dev_Display, LD * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!\n");
		mexPrintf("cudaMemcpy dev_Display failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	cudaFree(dev_R);
	cudaFree(dev_Display);


	mexPrintf("Exit FDK\n");
	return cudaStatus;
}