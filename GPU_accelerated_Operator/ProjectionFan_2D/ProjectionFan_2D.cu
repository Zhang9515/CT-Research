#include "ProjectionFan_2D.h"
#include <cstdlib>
// 2018/04/16
//__device__ const double PI = 3.141592653589793;
__device__ const double EPS = 1e-15;
__device__ const double ERR = 1e-5;

__global__ void ProjectionFan(const float *dev_Pic, double *dev_Projection, const double *dev_Pdomain, const double *dev_BetaScanRange,
	const double Center_t, const double Center_s, const double *dev_resolution, const int t_length, const int s_length, 
	const int Pstart, const int Betastart, const double Distance, const int LP, const int LBeta, const double *dev_RandomErr)
{  
	const unsigned short Pindex = threadIdx.x + Pstart;
	const unsigned short  Betaindex = blockIdx.x + Betastart;

	const unsigned long threadid = Betaindex * LP + Pindex;

	dev_Projection[threadid] = 0;

	double P = dev_Pdomain[Pindex];
	double Beta = dev_BetaScanRange[Betaindex];

	double resolution_1 = dev_resolution[0]; double resolution_2 = dev_resolution[1]; 

	// according to euler equation   
	double source_t = Center_t - Distance * sin(Beta);      // define the source in matlab coordinate
	//double source_t = Center_t - Distance * sin(Beta) + dev_RandomErr[threadid];  
	double source_s = Center_s + Distance * cos(Beta);

	double Theta = atan(P / Distance);        // radian angle in s'-t coordinate plane 

	double Smax = 2 * Distance;

	// define end detect point in matlab coordinate, Note that : 0 is the start
	double DetectPoint_tend = Center_t + Smax * sin(Theta) * cos(Beta) - (Distance - Smax * cos(Theta)) * sin(Beta);
	double DetectPoint_send = Center_s + Smax * sin(Theta) * sin(Beta) + (Distance - Smax * cos(Theta)) * cos(Beta);

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

	// to determine the range of t
	short t_signal = 0;

	if (DetectPoint_tend >= source_t)
		t_signal = 1;
	else
		t_signal = -1;

	// to determine the range of s
	short s_signal = 0;

	if (DetectPoint_send >= source_s)
		s_signal = 1;
	else
		s_signal = -1;

	// actual Size
	double tlow = 0, thigh = t_length*resolution_1, slow = 0, shigh = s_length*resolution_2;

	//compute the first and last point in the ROI
	// using DetectPoint_end set up projection equation
	double tlow_s = source_s + (tlow - source_t) * T2S;
	double thigh_s = source_s + (thigh - source_t) * T2S;

	double slow_t = source_t + (slow - source_s) * S2T;
	double shigh_t = source_t + (shigh - source_s) * S2T;

	//double *Range = new double [6];   //  XYXY small-big(number)
	double T1 = 0, S1 = 0, T2 = 0, S2 = 0;

	if (tlow_s >= 0 && tlow_s <= shigh )
	{
		T1 = tlow; S1 = tlow_s; 
		if (thigh_s >= 0 && thigh_s <= shigh)
		{
			T2 = thigh; S2 = thigh_s;
		}
		else if (slow_t >= 0 && slow_t <= thigh)
		{
			T2 = slow_t; S2 = slow;
		}
		else if (shigh_t >= 0 && shigh_t <= thigh)
		{
			T2 = shigh_t; S2 = shigh;
		}
	}
	else if (thigh_s >= 0 && thigh_s <= shigh )
	{
		T1 = thigh; S1 = thigh_s;
		if (slow_t >= 0 && slow_t <= thigh)
		{
			T2 = slow_t; S2 = slow;
		}
		else if (shigh_t >= 0 && shigh_t <= thigh)
		{
			T2 = shigh_t; S2 = shigh;
		}
	}
	else if (slow_t >= 0 && slow_t <= thigh)
	{
		T1 = slow_t; S1 = slow;
		if (shigh_t >= 0 && shigh_t <= thigh)
		{
			T2 = shigh_t; S2 = shigh; 
		}
	}
	else
	{
		//dev_Projection[threadid] = threadid;
		return;
	}

	// set the start point
	double TStart = 0, SStart = 0;
	if (Distancesq(T1, S1, source_t, source_s) >= Distancesq(T2, S2, source_t, source_s))
	{
		TStart = T2;
		SStart = S2;
	}
	else
	{
		TStart = T1;
		SStart = S1;
	}

	// adjust the order

	if (T2 < T1)
	{
		double c = T1;
		T1 = T2;
		T2 = c;
	}
	if (S2 < S1)
	{
		double c = S1;
		S1 = S2;
		S2 = c;
	}

	//// enter the ROI
	double weight = 0, Ray = 0;
	short GridT = 0, GridS = 0;        // candidate crosspoint index in matlab(0~t_length)
	double GridT_s = 0, GridS_t = 0;    // candidate crosspoint index in matlab(0~256)
	short DetectPoint_t = 0, DetectPoint_s = 0;   // current pixel index in matlab pixel index in matlab(0~255)
	long Pointid = 0;
	double TCross = TStart / resolution_1, SCross = SStart / resolution_2;     // current crosspoint index in matlab(0~256)
	
	//while (((XCross * dev_resolution[1]) >= Range[0]) && ((XCross * dev_resolution[1]) <= Range[2]) 
	//	&& ((YCross * dev_resolution[0]) >= Range[1]) && ((YCross * dev_resolution[0]) <= Range[3]))

	for (short i = 0;i<(t_length + s_length - 1);i++)
	{
		// judge whether XCross/YCross is integer
		if (abs(TCross - round(TCross)) < EPS)
		{
			GridT = round(TCross) + t_signal;
		}
		else
		{
			GridT = floor(TCross) + flag1to1or_1to0(t_signal);
		}
		GridT_s = (source_s + (GridT * resolution_1 - source_t) * T2S) / resolution_2;

		if (abs(SCross - round(SCross)) < EPS)
		{
			GridS = round(SCross) + s_signal;
		}
		else
		{
			GridS = floor(SCross) + flag1to1or_1to0(s_signal);
		}
		GridS_t = (source_t + (GridS * resolution_2 - source_s) * S2T) / resolution_1;

		//judge which crosspoint is the nearest, means the smallest distance
		if (Distancesq(GridT, GridT_s, TCross, SCross) <= Distancesq(GridS_t, GridS, TCross, SCross))
		{
			weight = sqrt(Distancesq(GridT * resolution_1, GridT_s * resolution_2, TCross * resolution_1, SCross * resolution_2));
			DetectPoint_t = floor(MID(GridT, TCross));                 // the midpoint locates the pixel
			DetectPoint_s = floor(MID(GridT_s, SCross));
			
			TCross = GridT;    // update
			SCross = GridT_s;		
		}
		else
		{
			weight = sqrt(Distancesq(GridS_t * resolution_1, GridS * resolution_2, TCross * resolution_1, SCross * resolution_2));
			DetectPoint_t = floor(MID(GridS_t, TCross));                 // the midpoint locates the pixel
			DetectPoint_s = floor(MID(GridS, SCross));
			TCross = GridS_t;    // update
			SCross = GridS;
		}

		//judge whether the point is in the ROI
		if ((DetectPoint_t >= 0) && (DetectPoint_t < t_length) && (DetectPoint_s >= 0) && (DetectPoint_s < s_length))
		{
			Pointid = DetectPoint_s * t_length + DetectPoint_t;
			Ray += weight * dev_Pic[Pointid];
		}
		else
		{
			//dev_Projection[threadid] = 9000;
			break;
		}

	}

	__syncthreads();
	dev_Projection[threadid] = Ray;
}

// Helser function for using CUDA to add vectors in parallel.
cudaError_t ProjectionFan_2D(const float *Pic, double *Projection, const double *BetaScanRange, const double *Pdomain,
	const int t_length, const int s_length, const double Center_t, const double Center_s, const int LBeta, const int LP,
	const double Distance, const double *resolution)
{
	mexPrintf("Hello GenMatParalell!\n");

	float *dev_Pic = 0;
	double *dev_Pdomain = 0, *dev_BetaScanRange = 0, *dev_Projection = 0, *dev_RandomErr = 0;
	double *dev_resolution = 0;

	int threadcubic_x = MIN(threadX, LP);
	int blockcubic_x = MIN(blockX, LBeta);
	int LPResidual = LP % threadX;
	int LBetaResidual = LBeta % blockX;
	int PTime = LP / threadX;
	int BetaTime = LBeta / blockX;
	int Pstart = 0;
	int Betastart = 0;
	double Beta = 0;

	const dim3 thread_cubic(threadcubic_x, 1, 1);
	const dim3 block_cubic(blockcubic_x, 1, 1);
	dim3 thread_cubic_residual(1, 1, 1);  // initial
	dim3 block_cubic_residual(1, 1, 1);  // initial

	mexPrintf("threadcubic_x: %d blockcubic_x: %d LPResidual: %d LBetaResidual: %d\n",
		threadcubic_x, blockcubic_x, LPResidual, LBetaResidual);

	if (LPResidual != 0)
	{
		thread_cubic_residual.x = LPResidual;
	}
	if (LBetaResidual != 0)
	{
		block_cubic_residual.x = LBetaResidual;
	}

	cudaError_t cudaStatus;

	double *RandomErr =new double[LBeta * LP];
	for (int beta = 0; beta < LBeta; beta++)
	{
		for (int P = 0; P < LP; P++)
		{
			RandomErr[beta * LP + P] = 0.01 * (rand()/double(RAND_MAX)*2-1);
		}
	}

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		mexPrintf("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
		goto Error;
	}

	mexPrintf("Call for GPU space.\n");

	// Allocate GPU buffers for three vectors (two input, one output).
	cudaStatus = cudaMalloc((void**)&dev_Pic, t_length * s_length * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pic cudaMalloc failed!");
		mexPrintf("dev_Pic cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Projection, LBeta * LP * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Projection cudaMalloc failed!");
		mexPrintf("dev_Projection cudaMalloc failed!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_RandomErr, LBeta * LP * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Projection cudaMalloc failed!");
		mexPrintf("dev_Projection cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Pdomain, LP * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_t_Range cudaMalloc failed!");
		mexPrintf("dev_t_Range cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_BetaScanRange, LBeta * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_thetaRange cudaMalloc failed!");
		mexPrintf("dev_thetaRange cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_resolution, 2 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution cudaMalloc failed!");
		mexPrintf("dev_resolution cudaMalloc failed!\n");
		goto Error;
	}

	mexPrintf("Copy data from CPU to GPU.\n");

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_Pic, Pic, t_length * s_length * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "thetaRange cudaMemcpy failed!");
		mexPrintf("thetaRange cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Pdomain, Pdomain, LP * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "t_Range cudaMemcpy failed!");
		mexPrintf("t_Range cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_BetaScanRange, BetaScanRange, LBeta * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "thetaRange cudaMemcpy failed!");
		mexPrintf("thetaRange cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_resolution, resolution, 2 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution cudaMemcpy failed!");
		mexPrintf("dev_resolution cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_RandomErr, RandomErr, LBeta * LP * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_RandomErr cudaMemcpy failed!");
		mexPrintf("dev_RandomErr cudaMemcpy failed!\n");
		goto Error;
	}

	mexPrintf("Launch computation projection of each lines.\n");
	
	// Launch a kernel on the GPU with one thread for each element.
	
	for (int numP = 0; numP < PTime; numP++)
	{
		for (int numB = 0; numB < BetaTime; numB++)
		{
			Pstart = numP * threadX;
			Betastart = numB * blockX;
			//mexPrintf("%d %d\n", Pstart, Betastart);
			ProjectionFan << <block_cubic, thread_cubic >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_BetaScanRange,
				Center_t, Center_s, dev_resolution, t_length, s_length, Pstart, Betastart, Distance, LP, LBeta, dev_RandomErr);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionFan launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionFan launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at Pstart: %d Betastart: %d \n", Pstart, Betastart);
				goto Error;
			}
		}		
	}
		
	if (LPResidual != 0)
	{
		Pstart = LP - LPResidual;			
		if (LBetaResidual != 0)
		{
			Betastart = LBeta - LBetaResidual;
			//("%d %d\n", Pstart, Betastart);
			ProjectionFan << <block_cubic_residual, thread_cubic_residual >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_BetaScanRange,
				Center_t, Center_s, dev_resolution, t_length, s_length, Pstart, Betastart, Distance, LP, LBeta, dev_RandomErr);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionFan launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionFan launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at Pstart: %d Betastart: %d \n", Pstart, Betastart);
				goto Error;
			}
		}

		for (int numB = 0; numB < BetaTime; numB++)
		{
			Betastart = numB * blockX;
			//("%d %d\n", Pstart, Betastart);
			ProjectionFan << <block_cubic, thread_cubic_residual >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_BetaScanRange,
				Center_t, Center_s, dev_resolution, t_length, s_length, Pstart, Betastart, Distance, LP, LBeta, dev_RandomErr);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionFan launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionFan launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at Pstart: %d Betastart: %d \n", Pstart, Betastart);
				goto Error;
			}
		}				
	}
	if (LBetaResidual != 0)
	{
		Betastart = LBeta - LBetaResidual;
		for (int numP = 0; numP < PTime; numP++)
		{
			Pstart = numP * threadX;
			//mexPrintf("%d %d\n", Pstart, Betastart);
			ProjectionFan << <block_cubic_residual, thread_cubic >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_BetaScanRange,
				Center_t, Center_s, dev_resolution, t_length, s_length, Pstart, Betastart, Distance, LP, LBeta, dev_RandomErr);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionFan launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionFan launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at Pstart: %d Betastart: %d \n", Pstart, Betastart);
				goto Error;
			}
		}		
	}


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Projection launch failed: %s\n", cudaGetErrorString(cudaStatus));
		mexPrintf("Projection launch failed\n");
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		mexPrintf("cudaDeviceSynchronize returned error %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(Projection, dev_Projection, LP * LBeta * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!\n");
		mexPrintf("cudaMemcpy failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	cudaFree(dev_Pdomain);
	cudaFree(dev_BetaScanRange);
	cudaFree(dev_Projection);
	cudaFree(dev_Pic);
	cudaFree(dev_resolution);
	cudaFree(dev_RandomErr);

	return cudaStatus;
}
