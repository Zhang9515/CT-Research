#include "ImplementUI.h"

#define pow2(x) (1.0*(x)*(x))
__device__ const double PI = 3.141592653589793;

__global__ void GFunction(double *dev_G, const double GamaInt, const double *dev_Xigamadomain, const double Distance_s2d)
{

	const unsigned int Gamaindex = threadIdx.x;
	const unsigned int Xigamaindex = blockIdx.x;
	double xigama = dev_Xigamadomain[Xigamaindex];

	double Proportion = Distance_s2d / sqrt(pow2(Distance_s2d) + pow2(xigama));
	double GamaIntP = GamaInt * Proportion;

	// S_L filter 
	if (Gamaindex == 0)
		dev_G[THREAD_SIZE_X - 1 + Xigamaindex * Filter_SIZE] = 2.0 / pow2(PI*GamaIntP);
	else
	{
		dev_G[(THREAD_SIZE_X - 1) - Gamaindex + Xigamaindex * Filter_SIZE] = -2.0 * pow2(Gamaindex) / pow2(PI) / (4 * pow2(Gamaindex) - 1) / pow2(sin(Gamaindex * GamaIntP));
		dev_G[(THREAD_SIZE_X - 1) + Gamaindex + Xigamaindex * Filter_SIZE] = -2.0 * pow2(Gamaindex) / pow2(PI) / (4 * pow2(Gamaindex) - 1) / pow2(sin(Gamaindex * GamaIntP));
	}
}

__global__ void PreWeightFiltration(double *dev_Rcov, double *dev_R, const double *dev_G, const double *dev_Gamadomain,
	const double *dev_BetaScanRange, const double *dev_Xigamadomain, const double betaStartAngle, const double MaxGama,
	const double Distance, const double Distance_s2d, const double GamaInt)
{
  
	//const unsigned int Gamaindex = blockIdx.x * blockDim.x + threadIdx.x;
	//const unsigned int Xigamaindex = blockIdx.y * blockDim.y + threadIdx.y;
	//const unsigned int Betaindex = blockIdx.z * blockDim.z + threadIdx.z;
	//const unsigned long thread_id = Betaindex * ( gridDim.x * gridDim.y * blockDim.x * blockDim.y ) 
	//	+ Xigamaindex * ( gridDim.x * blockDim.x ) + Gamaindex ;
	
	const unsigned int Gamaindex = threadIdx.x;
	const unsigned int Xigamaindex = blockIdx.x;
	const unsigned int Betaindex = blockIdx.y;
	const unsigned long thread_id = Betaindex * (gridDim.x * blockDim.x)
		+ Xigamaindex * blockDim.x + Gamaindex;

	double Weight = 0 ;

	double beta = dev_BetaScanRange[Betaindex];
	double gama = dev_Gamadomain[Gamaindex];
	double xigama = dev_Xigamadomain[Xigamaindex];

	double Proportion = Distance_s2d / sqrt(pow2(Distance_s2d) + pow2(xigama));
	double betaP = beta * Proportion;
	double gamaP = gama * Proportion;
	double GamaIntP = GamaInt * Proportion;
	double DistanceP = Distance / Proportion;
	double betaStartAngleP = betaStartAngle * Proportion;
	double MaxGamaP = MaxGama * Proportion;

	// Parker function
	if (betaP >= betaStartAngleP && betaP < betaStartAngleP + MaxGamaP - 2 * gamaP)
		Weight = pow2(sin(PI / 4 * betaP / (MaxGamaP / 2 - gamaP)));
	else if(betaP >= betaStartAngleP + MaxGamaP - 2 * gamaP && betaP <= betaStartAngleP + PI - 2 * gamaP)
		Weight = 1;
	else if(betaP > betaStartAngleP + PI - 2 * gamaP && betaP <= betaStartAngleP + PI + MaxGamaP)
		Weight = pow2(sin(PI / 4 * (PI + MaxGamaP - betaP) / (MaxGamaP / 2 + gamaP)));
	else
		Weight = 0;

	dev_R[thread_id] = dev_R[thread_id] * DistanceP * cos(gamaP) * Weight;    // directly cover the input

	__syncthreads();
	double Rcovsum = 0;
	// convolution
	for (int g=0 ; g<THREAD_SIZE_X; g++)
	{
		//convolution
		//dev_Rcov[thread_id] += dev_R[Betaindex * (gridDim.x * gridDim.y * blockDim.x * blockDim.y)
		//	+ Xigamaindex * (gridDim.x * blockDim.x) + g] * dev_G[(gridDim.x*blockDim.x - 1) + Gamaindex - g];
		Rcovsum += dev_R[Betaindex * (gridDim.x * blockDim.x)
			+ Xigamaindex * blockDim.x + g] * dev_G[(blockDim.x - 1) + Gamaindex - g + Xigamaindex * Filter_SIZE];
	}
	__syncthreads();

	dev_Rcov[thread_id] = GamaIntP * Rcovsum;

}

__global__ void BackProjection(const double *dev_Rcov, double *dev_Display, const int *dev_Size,
	const int z_length, const double Beta, const double Distance, const double Distance_s2d, 
	const double *dev_Gamadomain, const double *dev_Xigamadomain, const double GamaInt, const double XigamaInt, 
	const double BetaScanInt, const double minGama, const double maxGama, const double minXigama, const double maxXigama, 
	const int betaIndex, const int LGama, const int LXigama)
{
	const unsigned int Tindex = threadIdx.x;
	const unsigned int Sindex = blockIdx.x;
	const unsigned int Zindex = blockIdx.y;
	const unsigned long thread_id = Zindex * (gridDim.x * blockDim.x)
		+ Sindex * blockDim.x + Tindex;
	// initialize

	const double Resolution_z = 1.0 * dev_Size[2] / z_length;

	// define the source in ground coordinate
	double source_t = Distance * cos(Beta + PI / 2);            
	double source_s = Distance * sin(Beta + PI / 2);
	double source_z = 0;

	// rotation center
	double Center_t = dev_Size[0] / 2;                          
	double Center_s = dev_Size[1] / 2;
	double Center_z = dev_Size[2] / 2;

	double image_t = (Tindex - 0.5) - Center_t;  double image_s = (Sindex - 0.5) - Center_s; double image_z = (Zindex - 0.5) * Resolution_z - Center_z;           // image pixel in ground coordinate
	double L2 = pow2(image_t - source_t) + pow2(image_s - source_s) + pow2(image_z - source_z);
	
	// rotate in ground coordinate
	double dect_t = image_t * cos(Beta) + image_s * sin(Beta);
	double dect_s = -image_t * sin(Beta) + image_s * cos(Beta);
	double dect_z = image_z;
	
	// define the projection position on the detector
	double Xigama = Distance_s2d * tan(asin(dect_z / sqrt(L2)));
	double Gama = atan(dect_t / (Distance - dect_s));

	double Proportion = sqrt(L2) / sqrt(L2 - pow2(dect_z));

	unsigned int XigamaN1index = 0, XigamaN2index = 0, GamaN1index = 0, GamaN2index = 0;
	double Gama_domain1 = 0, Gama_domain2 = 0, Xigama_domain1 = 0, Xigama_domain2 = 0;
	double Xig1 = 0, Xig2 = 0, Gama1 = 0, Gama2 = 0;
	double Display_pBeta = 0;

	if ((Gama >= minGama) && (Gama < maxGama) && (Xigama >= minXigama) && (Xigama < maxXigama))
	{
		XigamaN1index = floor(fabs(Xigama - dev_Xigamadomain[0]) / XigamaInt);
		XigamaN2index = XigamaN1index + 1;
		GamaN1index = floor(fabs(Gama - dev_Gamadomain[0]) / GamaInt);
		GamaN2index = GamaN1index + 1;

		Gama_domain1 = dev_Gamadomain[GamaN1index]; Gama_domain2 = dev_Gamadomain[GamaN2index];
		Xigama_domain1 = dev_Xigamadomain[XigamaN1index]; Xigama_domain2 = dev_Xigamadomain[XigamaN2index];

		//bilinear interpolation
		Xig1 = fabs(Xigama - Xigama_domain1); Xig2 = fabs(Xigama_domain2 - Xigama);
		Gama1 = fabs(Gama - Gama_domain1); Gama2 = fabs(Gama_domain2 - Gama);

		//Display_pBeta = 1;
		Display_pBeta = (Xig2 * Gama2 * dev_Rcov[betaIndex * LGama * LXigama + XigamaN1index * LGama + GamaN1index]
			+ Xig1 * Gama2 * dev_Rcov[betaIndex * LGama * LXigama + XigamaN2index * LGama + GamaN1index] + Xig2 * Gama1 * dev_Rcov[betaIndex * LGama * LXigama + XigamaN1index * LGama + GamaN2index]
			+ Xig1 * Gama1 * dev_Rcov[betaIndex * LGama * LXigama + XigamaN2index * LGama + GamaN2index]) / (GamaInt * XigamaInt) / L2 *  (BetaScanInt / Proportion);

	}

	dev_Display[thread_id] += Display_pBeta;

}



// Helper function for using CUDA to add vectors in parallel.
cudaError_t FDKUI(double *Display, const double *R, const double *Xigamadomain, const double *Gamadomain,
	const double *BetaScanRange, const double betaStartAngle, const double MaxGama, const double Distance,
	const double Distance_s2d, const int LGama, const int LBeta, const int LXigama, const int *Size, const int z_length)
{
	mexPrintf("Hello FDKUI!\n");
	double *dev_R = 0, *dev_Rcov = 0;
	double *dev_G = 0;
	double *dev_Gamadomain = 0, *dev_Xigamadomain =0, *dev_BetaScanRange = 0;
	int *dev_Size = 0;

	double *dev_Display = 0;
	double GamaInt = fabs(Gamadomain[1] - Gamadomain[0]);
	double XigamaInt = fabs(Xigamadomain[1] - Xigamadomain[0]);
	double BetaScanInt = fabs(BetaScanRange[1] - BetaScanRange[0]);
	
	double maxGama = Gamadomain[0] > Gamadomain[LGama - 1] ? Gamadomain[0] : Gamadomain[LGama - 1];
	double minGama = Gamadomain[0] < Gamadomain[LGama - 1] ? Gamadomain[0] : Gamadomain[LGama - 1];
	double maxXigama = Xigamadomain[0] > Xigamadomain[LXigama - 1] ? Xigamadomain[0] : Xigamadomain[LXigama - 1];
	double minXigama = Xigamadomain[0] < Xigamadomain[LXigama - 1] ? Xigamadomain[0] : Xigamadomain[LXigama - 1];
	//mexPrintf("%lf %lf %lf %lf \n", maxGama, minGama, maxXigama, minXigama);

	const long LDisplay = Size[0] * Size[1] * z_length; 

	const dim3 thread_cubic(THREAD_SIZE_X, 1, 1);
	const dim3 block_cubic(BLOCK_SIZE_X, BLOCK_SIZE_Y, 1);
	
	const dim3 thread_cubic_Bp(Size[0], 1, 1);
	const dim3 block_cubic_Bp(Size[1], z_length, 1);

	cudaError_t cudaStatus;

	mexPrintf("start cuda\n");

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
		mexPrintf("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
        goto Error;
    }

	mexPrintf("call for space in GPU\n");

    // Allocate GPU buffers for three vectors (two input, one output).
    cudaStatus = cudaMalloc((void**)&dev_R, ARRAY_SIZE_IN_BYTES);
    if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_R cudaMalloc failed!\n");
		mexPrintf("dev_R cudaMalloc failed!\n");
        goto Error1;
    }

	cudaStatus = cudaMalloc((void**)&dev_Gamadomain, LGama * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Gamadomain cudaMalloc failed!\n");
		mexPrintf("dev_Gamadomain cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Xigamadomain, LXigama * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Xigamadomain cudaMalloc failed!\n");
		mexPrintf("dev_Xigamadomain cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_BetaScanRange, LBeta * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_BetaScanRange cudaMalloc failed!\n");
		mexPrintf("dev_BetaScanRange cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_G, LXigama * Filter_SIZE * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_G cudaMalloc failed!\n");
		mexPrintf("dev_G cudaMalloc failed!\n");
		goto Error1;
	}

	cudaStatus = cudaMalloc((void**)&dev_Rcov, ARRAY_SIZE_OUT_BYTES);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_BetaScanRange cudaMalloc failed!\n");
		mexPrintf("dev_BetaScanRange cudaMalloc failed!\n");
		goto Error;
	}

	//mexPrintf("copy data in CPU to GPU\n");

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_R, R, ARRAY_SIZE_IN_BYTES, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy R failed!\n");
		mexPrintf("cudaMemcpy R failed!\n");
        goto Error1;
    }

	cudaStatus = cudaMemcpy(dev_Gamadomain, Gamadomain, LGama * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Gamadomain failed!\n");
		mexPrintf("cudaMemcpy Gamadomain failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Xigamadomain, Xigamadomain, LXigama * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Xigamadomain failed!\n");
		mexPrintf("cudaMemcpy Xigamadomain failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_BetaScanRange, BetaScanRange, LBeta * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy BetaScanRange failed!\n");
		mexPrintf("cudaMemcpy v failed!\n");
		goto Error;
	}

	mexPrintf("start parallel computation\n");
	
	// Generate Filter
	GFunction << <LXigama, LGama >> >(dev_G, GamaInt, dev_Xigamadomain, Distance_s2d);

    // Preweight and filtration

	PreWeightFiltration << <block_cubic, thread_cubic >> >(dev_Rcov, dev_R, dev_G, dev_Gamadomain,
		dev_BetaScanRange, dev_Xigamadomain, betaStartAngle, MaxGama, Distance, Distance_s2d, GamaInt);

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
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error1;
	}

Error1:
	cudaFree(dev_R);
	cudaFree(dev_G);

	cudaStatus = cudaMalloc((void**)&dev_Display, LDisplay * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Display cudaMalloc failed!\n");
		mexPrintf("dev_Display cudaMalloc failed!\n");
		goto Error;
	}
	cudaMemset(dev_Display, 0, sizeof(double));

	cudaStatus = cudaMalloc((void**)&dev_Size, 3 * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Size cudaMalloc failed!\n");
		mexPrintf("dev_Size cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Size, Size, 3 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Size failed!\n");
		mexPrintf("cudaMemcpy Size failed!\n");
		goto Error;
	}

	//Backprojection
	for (int betaIndex = 0; betaIndex < LBeta; betaIndex++)
	{
		//int betaIndex = 0;
		BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_Rcov, dev_Display, dev_Size, z_length,
			BetaScanRange[betaIndex], Distance, Distance_s2d, dev_Gamadomain, dev_Xigamadomain, GamaInt, XigamaInt,
			BetaScanInt, minGama, maxGama, minXigama, maxXigama, betaIndex, LGama, LXigama);
	}
		
    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BackProjection launch failed: %s\n", cudaGetErrorString(cudaStatus));
		mexPrintf("BackProjection launch failed\n");
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(Display, dev_Display, LDisplay * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!\n");
		mexPrintf("cudaMemcpy dev_Display failed!\n");
        goto Error;
    }


Error:    
	cudaFree(dev_Gamadomain);
	cudaFree(dev_BetaScanRange);
	cudaFree(dev_Xigamadomain);
	cudaFree(dev_Rcov);
	cudaFree(dev_Display);
	cudaFree(dev_Size);

	mexPrintf("Exit FDKUI\n");
  return cudaStatus;
}
