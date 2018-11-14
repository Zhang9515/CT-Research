#include "ProjectionCone_3D.h"
#include <cstdlib>
// 2018/04/16
//__device__ const double PI = 3.141592653589793;
__device__ const double EPS = 1e-6;

__global__ void ProjectionCone(const float *dev_Pic, float *dev_Projection, const float *dev_Pdomain, const float *dev_Xigamadomain,
	const double Center_t, const double Center_s, const double Center_z, const double *dev_resolution,
	const int t_length, const int s_length, const int z_length, const float Beta, const int numB, const int Pstart,
	const int Xigamastart, const double Distance, const int LP, const int LXigama, const float *dev_RandomErr)
{  
	const unsigned int Pindex = threadIdx.x + Pstart;
	const unsigned int  Xigamaindex = blockIdx.x + Xigamastart;

	const unsigned long threadid = numB * LXigama * LP + Xigamaindex * LP + Pindex;

	dev_Projection[threadid] = 0;

	float P = dev_Pdomain[Pindex];
	float Xigama = dev_Xigamadomain[Xigamaindex];

	// according to euler equation   
	double source_t = Center_t - Distance * sin(Beta);      // define the source in matlab coordinate
	//double source_t = Center_t - Distance * sin(Beta) + dev_RandomErr[threadid];  
	double source_s = Center_s + Distance * cos(Beta);
	double source_z = Center_z;
	//double source_z = Center_z + dev_RandomErr[threadid];

	double Gama = atan(Xigama / Distance);          // radian angle in s - z coordinate plane
	double Distance_shift = Distance / cos(Gama);    // length of DO'

	double Theta = atan(P / Distance_shift);        // radian angle in s'-t coordinate plane 

	double Smax = 2 * Distance;

	// define end detect point in matlab coordinate, Note that : 0 is the start
	double DetectPoint_tend = Center_t + Smax * sin(Theta) * cos(Beta) - (Distance - Smax * cos(Theta) * cos(Gama)) * sin(Beta);
	double DetectPoint_send = Center_s + Smax * sin(Theta) * sin(Beta) + (Distance - Smax * cos(Theta) * cos(Gama)) * cos(Beta);
	double DetectPoint_zend = Center_z + Smax * cos(Theta) * sin(Gama);

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

	// to determine the range of z
	short z_signal = 0;

	if (DetectPoint_zend >= source_z)
		z_signal = 1;
	else
		z_signal = -1;

	// actual Size
	double tlow = 0, thigh = t_length*dev_resolution[0], slow = 0, shigh = s_length*dev_resolution[1],
		zlow = 0, zhigh = z_length*dev_resolution[2];

	//compute the first and last point in the ROI
	// using DetectPoint_end set up projection equation
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

	//double *Range = new double [6];   //  XYXY small-big(number)
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
		else if (zlow_t >= 0 && zlow_t <= thigh && zlow_s>= 0 && zlow_s <= shigh)
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

	// set the start point
	double TStart = 0, SStart = 0, ZStart = 0;
	if (Distancesq(T1, S1, Z1, source_t, source_s, source_z) >= Distancesq(T2, S2, Z2, source_t, source_s, source_z))
	{
		TStart = T2;
		SStart = S2;
		ZStart = Z2;
	}
	else
	{
		TStart = T1;
		SStart = S1;
		ZStart = Z1;
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
	if (Z2 < Z1)
	{
		double c = Z1;
		Z1 = Z2;
		Z1 = c;
	}

	//// enter the ROI
	double weight = 0, Ray = 0;
	int GridT = 0, GridS = 0, GridZ = 0;        // candidate crosspoint index in matlab(0~t_length)
	double GridT_s = 0, GridT_z = 0,
		GridS_t = 0, GridS_z = 0,
		GridZ_t = 0, GridZ_s = 0;    // candidate crosspoint index in matlab(0~256)
	int DetectPoint_t = 0, DetectPoint_s = 0, DetectPoint_z = 0;   // current pixel index in matlab pixel index in matlab(0~255)
	long Pointid = 0;
	double TCross = TStart / dev_resolution[0], SCross = SStart / dev_resolution[1],
		ZCross = ZStart / dev_resolution[2];     // current crosspoint index in matlab(0~256)
	
	int i = 0;
	//while (((XCross * dev_resolution[1]) >= Range[0]) && ((XCross * dev_resolution[1]) <= Range[2]) 
	//	&& ((YCross * dev_resolution[0]) >= Range[1]) && ((YCross * dev_resolution[0]) <= Range[3]))
	while (i < (t_length + s_length + z_length -2))
	{
		i++;
		// judge whether XCross/YCross is integer
		if (TCross - (double)((int)TCross) < EPS)
		{
			GridT = TCross + t_signal;
		}
		else
		{
			GridT = floor(TCross) + flag1to1or_1to0(t_signal);
		}
		GridT_s = (source_s + (GridT * dev_resolution[0] - source_t) * (DetectPoint_send - source_s) / (DetectPoint_tend - source_t)) / dev_resolution[1];
		GridT_z = (source_z + (GridT * dev_resolution[0] - source_t) * (DetectPoint_zend - source_z) / (DetectPoint_tend - source_t)) / dev_resolution[2];

		if (SCross - (double)((int)SCross) < EPS)
		{
			GridS = SCross + s_signal;
		}
		else
		{
			GridS = floor(SCross) + flag1to1or_1to0(s_signal);
		}
		GridS_t = (source_t + (GridS * dev_resolution[1] - source_s) * (DetectPoint_tend - source_t) / (DetectPoint_send - source_s)) / dev_resolution[0];
		GridS_z = (source_z + (GridS * dev_resolution[1] - source_s) * (DetectPoint_zend - source_z) / (DetectPoint_send - source_s)) / dev_resolution[2];

		if (ZCross - (double)((int)ZCross) < EPS)
		{
			GridZ = ZCross + z_signal;
		}
		else
		{
			GridZ = floor(ZCross) + flag1to1or_1to0(z_signal);
		}
		GridZ_t = (source_t + (GridZ * dev_resolution[2] - source_z) * (DetectPoint_tend - source_t) / (DetectPoint_zend - source_z)) / dev_resolution[0];
		GridZ_s = (source_s + (GridZ * dev_resolution[2] - source_z) * (DetectPoint_send - source_s) / (DetectPoint_zend - source_z)) / dev_resolution[1];

		//judge which crosspoint is the nearest, means the smallest distance
		if (Distancesq(GridT, GridT_s, GridT_z, TCross, SCross, ZCross) <= Distancesq(GridS_t, GridS, GridS_z, TCross, SCross, ZCross))
		{
			if (Distancesq(GridZ_t, GridZ_s, GridZ, TCross, SCross, ZCross) <= Distancesq(GridT, GridT_s, GridT_z, TCross, SCross, ZCross))
			{
				weight = sqrt(Distancesq(GridZ_t * dev_resolution[0], GridZ_s * dev_resolution[1],
					GridZ * dev_resolution[2], TCross * dev_resolution[0], SCross * dev_resolution[1], 
					ZCross * dev_resolution[2]));
				DetectPoint_t = floor(MID(GridZ_t, TCross));                 // the midpoint locates the pixel
				DetectPoint_s = floor(MID(GridZ_s, SCross));
				DetectPoint_z = floor(MID(GridZ, ZCross));
				TCross = GridZ_t;    // update
				SCross = GridZ_s;
				ZCross = GridZ;
			}
			else
			{
				weight = sqrt(Distancesq(GridT * dev_resolution[0], GridT_s * dev_resolution[1],
					GridT_z * dev_resolution[2], TCross * dev_resolution[0], SCross * dev_resolution[1],
					ZCross * dev_resolution[2]));
				DetectPoint_t = floor(MID(GridT, TCross));                 // the midpoint locates the pixel
				DetectPoint_s = floor(MID(GridT_s, SCross));
				DetectPoint_z = floor(MID(GridT_z, ZCross));
				TCross = GridT;    // update
				SCross = GridT_s;
				ZCross = GridT_z;
			}			
		}
		else
		{
			if (Distancesq(GridZ_t, GridZ_s, GridZ, TCross, SCross, ZCross) <= Distancesq(GridS_t, GridS, GridS_z, TCross, SCross, ZCross))
			{
				weight = sqrt(Distancesq(GridZ_t * dev_resolution[0], GridZ_s * dev_resolution[1],
					GridZ * dev_resolution[2], TCross * dev_resolution[0], SCross * dev_resolution[1],
					ZCross * dev_resolution[2]));
				DetectPoint_t = floor(MID(GridZ_t, TCross));                 // the midpoint locates the pixel
				DetectPoint_s = floor(MID(GridZ_s, SCross));
				DetectPoint_z = floor(MID(GridZ, ZCross));
				TCross = GridZ_t;    // update
				SCross = GridZ_s;
				ZCross = GridZ;
			}
			else
			{
				weight = sqrt(Distancesq(GridS_t * dev_resolution[0], GridS * dev_resolution[1],
					GridS_z * dev_resolution[2], TCross * dev_resolution[0], SCross * dev_resolution[1],
					ZCross * dev_resolution[2]));
				DetectPoint_t = floor(MID(GridS_t, TCross));                 // the midpoint locates the pixel
				DetectPoint_s = floor(MID(GridS, SCross));
				DetectPoint_z = floor(MID(GridS_z, ZCross));
				TCross = GridS_t;    // update
				SCross = GridS;
				ZCross = GridS_z;
			}
		}

		//judge whether the point is in the ROI
		if ((DetectPoint_t >= 0) && (DetectPoint_t <= (t_length - 1)) && (DetectPoint_s >= 0) && (DetectPoint_s <= (s_length - 1))
			&& (DetectPoint_z >= 0) && (DetectPoint_z <= (z_length - 1)))
		{
			Pointid = DetectPoint_z * t_length * s_length + DetectPoint_s * t_length + DetectPoint_t;
			Ray += weight * dev_Pic[Pointid];
		}
		else
		{
			//dev_Projection[threadid] = 9000;
			break;
		}

	}

	//__syncthreads();
	dev_Projection[threadid] = Ray;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t ProjectionCone_3D(const float *Pic, float *Projection, const float *BetaScanRange, const float *Pdomain,
	const float *Xigamadomain, const int t_length, const int s_length, const int z_length, const double Center_t,
	const double Center_s, const double Center_z, const int LBeta, const int LP, const int LXigama, const double Distance,
	const double *resolution)
{
	mexPrintf("Hello GenMatParalell!\n");

	float *dev_Pic = 0, *dev_Pdomain = 0, *dev_Xigamadomain = 0, *dev_Projection = 0, *dev_RandomErr = 0;
	double *dev_resolution = 0;

	int threadcubic_x = MIN(threadX, LP);
	int blockcubic_x = MIN(blockX, LXigama);
	int LPResidual = LP % threadX;
	int LXigamaResidual = LXigama % blockX;
	int PTime = LP / threadX;
	int XigamaTime = LXigama / blockX;
	int Pstart = 0;
	int Xigamastart = 0;
	float Beta = 0;
	const dim3 thread_cubic(threadcubic_x, 1, 1);
	const dim3 block_cubic(blockcubic_x, 1, 1);
	dim3 thread_cubic_residual(1, 1, 1);  // initial
	dim3 block_cubic_residual(1, 1, 1);  // initial

	mexPrintf("%d %d %d %d\n", threadcubic_x, blockcubic_x, LPResidual, LXigamaResidual);

	if (LPResidual != 0)
	{
		thread_cubic_residual.x = LPResidual;
	}
	if (LXigamaResidual != 0)
	{
		block_cubic_residual.x = LXigamaResidual;
	}

	cudaError_t cudaStatus;

	float *RandomErr =new float[LBeta * LP * LXigama];
	for (int beta = 0; beta < LBeta; beta++)
	{
		for (int P = 0; P < LP; P++)
		{
			for (int Xigama = 0; Xigama < LXigama; Xigama++)
			{
				RandomErr[beta * LXigama * LP + Xigama * LP + P] = 0.01 * (rand()/float(RAND_MAX)*2-1);
			}
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
	cudaStatus = cudaMalloc((void**)&dev_Pic, t_length * s_length * z_length * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pic cudaMalloc failed!");
		mexPrintf("dev_Pic cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Projection, LBeta * LP * LXigama * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Projection cudaMalloc failed!");
		mexPrintf("dev_Projection cudaMalloc failed!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_RandomErr, LBeta * LP * LXigama * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Projection cudaMalloc failed!");
		mexPrintf("dev_Projection cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Pdomain, LP * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_t_Range cudaMalloc failed!");
		mexPrintf("dev_t_Range cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Xigamadomain, LXigama * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution cudaMalloc failed!");
		mexPrintf("dev_resolution cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_resolution, 3 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution cudaMalloc failed!");
		mexPrintf("dev_resolution cudaMalloc failed!\n");
		goto Error;
	}

	mexPrintf("Copy data from CPU to GPU.\n");

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_Pic, Pic, t_length * s_length * z_length * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "thetaRange cudaMemcpy failed!");
		mexPrintf("thetaRange cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Pdomain, Pdomain, LP * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "t_Range cudaMemcpy failed!");
		mexPrintf("t_Range cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_Xigamadomain, Xigamadomain, LXigama * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "t_Range cudaMemcpy failed!");
		mexPrintf("t_Range cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_resolution, resolution, 3 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution cudaMemcpy failed!");
		mexPrintf("dev_resolution cudaMemcpy failed!\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_RandomErr, RandomErr, LBeta * LP * LXigama * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_RandomErr cudaMemcpy failed!");
		mexPrintf("dev_RandomErr cudaMemcpy failed!\n");
		goto Error;
	}
	dev
	mexPrintf("Launch computation projection of each lines.\n");
	
	// Launch a kernel on the GPU with one thread for each element.
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
				ProjectionCone << <block_cubic, thread_cubic >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_Xigamadomain,
					Center_t, Center_s, Center_z, dev_resolution, t_length, s_length, z_length, Beta, numB, Pstart, Xigamastart,
					Distance, LP, LXigama, dev_RandomErr);
			}		
		}
		
		if (LPResidual != 0)
		{
			Pstart = LP - LPResidual;			
			if (LXigamaResidual != 0)
			{
				Xigamastart = LXigama - LXigamaResidual;
				//("%d %d\n", Pstart, Xigamastart);
				ProjectionCone << <block_cubic_residual, thread_cubic_residual >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_Xigamadomain,
					Center_t, Center_s, Center_z, dev_resolution, t_length, s_length, z_length, Beta, numB, Pstart, Xigamastart,
					Distance, LP, LXigama, dev_RandomErr);
			}

			for (int numX = 0; numX < XigamaTime; numX++)
			{
				Xigamastart = numX * blockX;
				//("%d %d\n", Pstart, Xigamastart);
				ProjectionCone << <block_cubic, thread_cubic_residual >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_Xigamadomain,
					Center_t, Center_s, Center_z, dev_resolution, t_length, s_length, z_length, Beta, numB, Pstart, Xigamastart,
					Distance, LP, LXigama, dev_RandomErr);
			}				
		}
		if (LXigamaResidual != 0)
		{
			Xigamastart = LXigama - LXigamaResidual;
			for (int numP = 0; numP < PTime; numP++)
			{
				Pstart = numP * threadX;
				//mexPrintf("%d %d\n", Pstart, Xigamastart);
				ProjectionCone << <block_cubic_residual, thread_cubic >> >(dev_Pic, dev_Projection, dev_Pdomain, dev_Xigamadomain,
					Center_t, Center_s, Center_z, dev_resolution, t_length, s_length, z_length, Beta, numB, Pstart, Xigamastart,
					Distance, LP, LXigama, dev_RandomErr);
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
	cudaStatus = cudaMemcpy(Projection, dev_Projection, LBeta * LP * LXigama * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!\n");
		mexPrintf("cudaMemcpy failed! %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	cudaFree(dev_Pdomain);
	cudaFree(dev_Xigamadomain);
	cudaFree(dev_Projection);
	cudaFree(dev_Pic);
	cudaFree(dev_resolution);

	return cudaStatus;
}
