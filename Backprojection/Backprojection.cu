#include "Backprojection.h"
// 2018/11/16 apply GPU acceleration

__device__ const double PI = 3.141592653589793;

// thiss
__global__ void BackProjection(const float *dev_R, float *dev_Display, const double *dev_Size,
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
	// image pixel in ground coordinate
	double image_t = (Tindex + 0.5) * Resolution_t - Center_t;  
	double image_s = (Sindex + 0.5) * Resolution_s - Center_s; 
	double image_z = (Zindex + 0.5) * Resolution_z - Center_z;

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
	double backweight = 0;
	double LengthinROI = 0;

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

		Display_pBeta = (Xig2 * P2 * dev_R[betaIndex * LP * LXigama + XigamaN1index * LP + PN1index]
			+ Xig1 * P2 * dev_R[betaIndex * LP * LXigama + XigamaN2index * LP + PN1index] + Xig2 * P1 * dev_R[betaIndex * LP * LXigama + XigamaN1index * LP + PN2index]
			+ Xig1 * P1 * dev_R[betaIndex * LP * LXigama + XigamaN2index * LP  + PN2index]) / (PInt * XigamaInt);

	//	 the way to compute backweight is to get the cross length in the specific pixel and the whole ROI

	//	 according to euler equation   
		double source_t = Center_t - Distance * sin(Beta);      // define the source in matlab coordinate														
		double source_s = Center_s + Distance * cos(Beta);
		double source_z = Center_z;
		
		// in matlab coordinate
		//  assume the projection line go through the center of the current pixel 
		double DetectPoint_tend = image_t + Center_t;
		double DetectPoint_send = image_s + Center_s;
		double DetectPoint_zend = image_z + Center_z;

	//	first compute length in whole ROI
	//	 actual Size
		double tlow = 0, thigh = t_length * Resolution_t, slow = 0, shigh = s_length * Resolution_s,
			zlow = 0, zhigh = z_length * Resolution_z;

	//	compute the first and last point in the ROI
	//	 using DetectPoint_end set up projection equation
		double tlow_s = source_s + (tlow - source_t) * (DetectPoint_send - source_s) / (DetectPoint_tend - source_t);
		double tlow_z = source_z + (tlow - source_t) * (DetectPoint_zend - source_z) / (DetectPoint_tend - source_t);
		double thigh_s = source_s + (thigh - source_t) * (DetectPoint_send - source_s) / (DetectPoint_tend - source_t);
		double thigh_z = source_z + (thigh - source_t) * (DetectPoint_zend - source_z) / (DetectPoint_tend - source_t);

		double slow_t = source_t + (slow - source_s) * (DetectPoint_tend - source_t) / (DetectPoint_send - source_s);
		double slow_z = source_z + (slow - source_s) * (DetectPoint_zend - source_z) / (DetectPoint_send - source_s);
		double shigh_t = source_t + (shigh - source_s) * (DetectPoint_tend - source_t) / (DetectPoint_send - source_s);
		double shigh_z = source_z + (shigh - source_s) * (DetectPoint_zend - source_z) / (DetectPoint_send - source_s);

		double zlow_t = source_t + (zlow - source_z) * (DetectPoint_tend - source_t) / (DetectPoint_zend - source_z);
		double zlow_s = source_s + (zlow - source_z) * (DetectPoint_send - source_s) / (DetectPoint_zend - source_z);
		double zhigh_t = source_t + (zhigh - source_z) * (DetectPoint_tend - source_t) / (DetectPoint_zend - source_z);
		double zhigh_s = source_s + (zhigh - source_z) * (DetectPoint_send - source_s) / (DetectPoint_zend - source_z);

	//	double *Range = new double [6];   //  XYXY small-big(number)
		double T1 = 0, S1 = 0, Z1 = 0, T2 = 0, S2 = 0, Z2 = 0;

		if (tlow_s >= 0 && tlow_s <= shigh && tlow_z >= 0 && tlow_z <= zhigh)
		{
			T1 = tlow; S1 = tlow_s; Z1 = tlow_z;
			if (thigh_s >= 0 && thigh_s <= shigh && thigh_z >= 0 && thigh_z <= zhigh)
			{
				T2 = thigh; S2 = thigh_s; Z2 = thigh_z;
			}
			else if (slow_t >= 0 && slow_t <= thigh && slow_z >= 0 && slow_z <= zhigh)
			{
				T2 = slow_t; S2 = slow; Z2 = slow_z;
			}
			else if (shigh_t >= 0 && shigh_t <= thigh && shigh_z >= 0 && shigh_z <= zhigh)
			{
				T2 = shigh_t; S2 = shigh; Z2 = shigh_z;
			}
			else if (zlow_t >= 0 && zlow_t <= thigh && zlow_s >= 0 && zlow_s <= shigh)
			{
				T2 = zlow_t; S2 = zlow_s; Z2 = zlow;
			}
			else if (zhigh_t >= 0 && zhigh_t <= thigh && zhigh_s >= 0 && zhigh_s <= shigh)
			{
				T2 = zhigh_t; S2 = zhigh_s; Z2 = zhigh;
			}
		}
		else if (thigh_s >= 0 && thigh_s <= shigh && thigh_z >= 0 && thigh_z <= zhigh)
		{
			T1 = thigh; S1 = thigh_s; Z1 = thigh_z;
			if (slow_t >= 0 && slow_t <= thigh && slow_z >= 0 && slow_z <= zhigh)
			{
				T2 = slow_t; S2 = slow; Z2 = slow_z;
			}
			else if (shigh_t >= 0 && shigh_t <= thigh && shigh_z >= 0 && shigh_z <= zhigh)
			{
				T2 = shigh_t; S2 = shigh; Z2 = shigh_z;
			}
			else if (zlow_t >= 0 && zlow_t <= thigh && zlow_s >= 0 && zlow_s <= shigh)
			{
				T2 = zlow_t; S2 = zlow_s; Z2 = zlow;
			}
			else if (zhigh_t >= 0 && zhigh_t <= thigh && zhigh_s >= 0 && zhigh_s <= shigh)
			{
				T2 = zhigh_t; S2 = zhigh_s; Z2 = zhigh;
			}
		}
		else if (slow_t >= 0 && slow_t <= thigh && slow_z >= 0 && slow_z <= zhigh)
		{
			T1 = slow_t; S1 = slow; Z1 = slow_z;
			if (shigh_t >= 0 && shigh_t <= thigh && shigh_z >= 0 && shigh_z <= zhigh)
			{
				T2 = shigh_t; S2 = shigh; Z2 = shigh_z;
			}
			else if (zlow_t >= 0 && zlow_t <= thigh && zlow_s >= 0 && zlow_s <= shigh)
			{
				T2 = zlow_t; S2 = zlow_s; Z2 = zlow;
			}
			else if (zhigh_t >= 0 && zhigh_t <= thigh && zhigh_s >= 0 && zhigh_s <= shigh)
			{
				T2 = zhigh_t; S2 = zhigh_s; Z2 = zhigh;
			}
		}
		else if (shigh_t >= 0 && shigh_t <= thigh && shigh_z >= 0 && shigh_z <= zhigh)
		{
			T1 = shigh_t; S1 = shigh; Z1 = shigh_z;
			if (zlow_t >= 0 && zlow_t <= thigh && zlow_s >= 0 && zlow_s <= shigh)
			{
				T2 = zlow_t; S2 = zlow_s; Z2 = zlow;
			}
			else if (zhigh_t >= 0 && zhigh_t <= thigh && zhigh_s >= 0 && zhigh_s <= shigh)
			{
				T2 = zhigh_t; S2 = zhigh_s; Z2 = zhigh;
			}
		}
		else if (zlow_t >= 0 && zlow_t <= thigh && zlow_s >= 0 && zlow_s <= shigh)
		{
			T1 = zlow_t; S1 = zlow_s; Z1 = zlow;
			if (zhigh_t >= 0 && zhigh_t <= thigh && zhigh_s >= 0 && zhigh_s <= shigh)
			{
				T2 = zhigh_t; S2 = zhigh_s; Z2 = zhigh;
			}
		}
		else
		{
			//dev_Projection[threadid] = threadid;
			return;
		}

		LengthinROI = Distance(T1, S1, Z1, T2, S2, Z2);

	//	secondly compute length in a single pixel, the process is very similar to the previous.
	//	because this time the line goes through the center point of the pixel. So some kind of symmetry happens.
	//	since the global function can not call exterior function, so the previous code will be copied here.

	//	 actual Size in matlab coordinate
		tlow = Tindex * Resolution_t; thigh = (Tindex + 1) * Resolution_t; 
		slow = Sindex * Resolution_s; shigh = (Sindex + 1) * Resolution_s;
		zlow = Zindex * Resolution_z; zhigh = (Zindex + 1) * Resolution_z;

	//	compute the first and last point in the ROI
	//	 using DetectPoint_end set up projection equation

		tlow_s = source_s + (tlow - source_t) * (DetectPoint_send - source_s) / (DetectPoint_tend - source_t);
		tlow_z = source_z + (tlow - source_t) * (DetectPoint_zend - source_z) / (DetectPoint_tend - source_t);

		slow_t = source_t + (slow - source_s) * (DetectPoint_tend - source_t) / (DetectPoint_send - source_s);
		slow_z = source_z + (slow - source_s) * (DetectPoint_zend - source_z) / (DetectPoint_send - source_s);

		zlow_t = source_t + (zlow - source_z) * (DetectPoint_tend - source_t) / (DetectPoint_zend - source_z);
		zlow_s = source_s + (zlow - source_z) * (DetectPoint_send - source_s) / (DetectPoint_zend - source_z);

	//	double *Range = new double [6];   //  XYXY small-big(number)
		T1 = 0; S1 = 0; Z1 = 0; T2 = 0; S2 = 0; Z2 = 0;

		if (tlow_s >= slow && tlow_s <= shigh && tlow_z >= zlow && tlow_z <= zhigh)
		{
			T1 = tlow; S1 = tlow_s; Z1 = tlow_z;
			T2 = 2 * DetectPoint_tend - T1; 
			S2 = 2 * DetectPoint_send - S1; 
			Z2 = 2 * DetectPoint_zend - Z1;
	//		symmetry
		}
		else if (slow_t >= tlow && slow_t <= thigh && slow_z >= zlow && slow_z <= zhigh)
		{
			T1 = slow_t; S1 = slow; Z1 = slow_z;
			T2 = 2 * DetectPoint_tend - T1; 
			S2 = 2 * DetectPoint_send - S1; 
			Z2 = 2 * DetectPoint_zend - Z1;
	//		symmetry
		}
		else if (zlow_t >= tlow && zlow_t <= thigh && zlow_s >= slow && zlow_s <= shigh)
		{
			T1 = zlow_t; S1 = zlow_s; Z1 = zlow;
			T2 = 2 * DetectPoint_tend - T1; 
			S2 = 2 * DetectPoint_send - S1; 
			Z2 = 2 * DetectPoint_zend - Z1;
	//		symmetry		
		}
		else
		{
	//		dev_Projection[threadid] = threadid;
			return;
		}

		double LengthinPixel = Distance(T1, S1, Z1, T2, S2, Z2);
		backweight = LengthinPixel / LengthinROI ; 
	}

	dev_Display[thread_id] += Display_pBeta * backweight ;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t BackPro(float *Display, const float *R, const float *Xigamadomain, const float *Pdomain,
	const float *BetaScanRange, const double Distance, const int LBeta, const int LP, const int LXigama,
	const double *Size, const int t_length, const int s_length, const int z_length)
{
	mexPrintf("BackPro!\n");
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

	// Allocate GPU buffers for four vectors (4 inputs).

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

	cudaStatus = cudaMalloc((void**)&dev_BetaScanRange, LBeta * sizeof(float));
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
		goto Error;
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
		mexPrintf("cudaMemcpy BetaScanRange failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Xigamadomain, Xigamadomain, LXigama * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy Xigamadomain failed!\n");
		mexPrintf("cudaMemcpy Xigamadomain failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	mexPrintf("start parallel computation\n");
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
	
    // output
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
		BackProjection << <block_cubic_Bp, thread_cubic_Bp >> > (dev_R, dev_Display, dev_Size, t_length, s_length, z_length,
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
	cudaFree(dev_R);
	cudaFree(dev_BetaScanRange);
	cudaFree(dev_Pdomain);
	cudaFree(dev_Xigamadomain);
	cudaFree(dev_Display);
	cudaFree(dev_Size);

	mexPrintf("Exit Bakprojection\n");
	return cudaStatus;
}