#include "ProjectionParallel_2D.h"
// 2018/04/09
// 2019/03/24 edited by ZXZ
__device__ const double PI = 3.141592653589793;
__device__ const double EPS = 1e-10;

__global__ void ProjectionParallel(const float *dev_Pic, double *dev_Projection, const float *dev_t_Range, const float *dev_thetaRange,
	const double Center_y, const double Center_x, const double *dev_resolution, const int height, const int width, const int Lt, 
	const int tstart, const int thetastart)
{	
	const short tindex = threadIdx.x + tstart;
	const short thetaindex = blockIdx.x + thetastart;
	//for debug{
	//const int tindex = 364;
	//const int thetaindex = 2343;
	//}
	const int threadid = thetaindex * Lt + tindex;

	dev_Projection[threadid] = 0;

	float t = dev_t_Range[tindex];
	float theta = dev_thetaRange[thetaindex];
	double Smax = MAX(fabs(dev_t_Range[0]), fabs(dev_t_Range[Lt-1]));

	// according to euler equation
	double DetectPoint_xstart = Center_x + t * cos(-theta) + (-Smax) * sin(-theta);
	double DetectPoint_ystart = Center_y - t * sin(-theta) + (-Smax) * cos(-theta);
	// define end detect point in matlab coordinate, Note that : 0 is the start
	double DetectPoint_xend = Center_x + t * cos(-theta) + Smax * sin(-theta);
	double DetectPoint_yend = Center_y - t * sin(-theta) + Smax * cos(-theta);
	
	double X2Y = (DetectPoint_yend - DetectPoint_ystart) / (DetectPoint_xend - DetectPoint_xstart + EPS);
	double Y2X = 1 / (X2Y + EPS);

	// to determine the range of y
	short y_signal = 0;

	if (DetectPoint_yend >= DetectPoint_ystart)
		y_signal = 1;
	else if (DetectPoint_yend < DetectPoint_ystart)
		y_signal = -1;

	// to determine the range of x
	short x_signal = 0;

	if (DetectPoint_xend >= DetectPoint_xstart)			
		x_signal = 1;			
	else if (DetectPoint_xend < DetectPoint_xstart)
		x_signal = -1;			

	
	// // actual Size
	double Xlow = 0, Xhigh = width*dev_resolution[1], Ylow = 0, Yhigh = height*dev_resolution[0];

	//compute the first and last point in the ROI
	double Ylow_x = DetectPoint_xstart + (Ylow - DetectPoint_ystart) * Y2X /*/ (tan(theta + PI / 2) + EPS)*/;
	double Yhigh_x = DetectPoint_xstart + (Yhigh - DetectPoint_ystart) * Y2X /*/ (tan(theta + PI / 2) + EPS)*/;
	
	double Xlow_y = DetectPoint_ystart + (Xlow - DetectPoint_xstart) * X2Y /** tan(theta + PI / 2)*/;
	double Xhigh_y = DetectPoint_ystart + (Xhigh - DetectPoint_xstart) * X2Y /** tan(theta + PI / 2)*/;
	//double *Range = new double [4];   //  XYXY small-big(number)
	double X1 = 0, Y1 = 0, X2 = 0, Y2 = 0;

	if (Ylow_x >= 0 && Ylow_x <= Xhigh)
	{
		X1 = Ylow_x; Y1 = Ylow;
		if (Xlow_y >= 0 && Xlow_y <= Yhigh)
		{
			X2 = Xlow; Y2 = Xlow_y;
		}

		else if (Xhigh_y >= 0 && Xhigh_y <= Yhigh)
		{
			X2 = Xhigh; Y2 = Xhigh_y;
		}
		else if (Yhigh_x >= 0 && Yhigh_x <= Xhigh)
		{
			X2 = Yhigh_x; Y2 = Yhigh;
		}

	}
	else if (Yhigh_x >= 0 && Yhigh_x <= Xhigh)
	{
		X1 = Yhigh_x; Y1 = Yhigh;
		if (Xlow_y >= 0 && Xlow_y <= Yhigh)
		{
			X2 = Xlow; Y2 = Xlow_y;
		}
		else if (Xhigh_y >= 0 && Xhigh_y <= Yhigh)
		{
			X2 = Xhigh; Y2 = Xhigh_y;
		}

	}
	else if (Xlow_y >= 0 && Xlow_y <= Yhigh)
	{
		X1 = Xlow; Y1 = Xlow_y;
		if (Xhigh_y >= 0 && Xhigh_y <= Yhigh)
		{
			X2 = Xhigh; Y2 = Xhigh_y;
		}

	}
	else
	{
		//dev_Projection[threadid] = 8000;
		return;
	}
		
	// set the start point
	double XStart = 0, YStart = 0;
	if (Distancesq(X1, Y1, DetectPoint_xstart, DetectPoint_ystart) >= Distancesq(X2, Y2, DetectPoint_xstart, DetectPoint_ystart))
	{
		XStart = X2;
		YStart = Y2;
	}
	else
	{
		XStart = X1;
		YStart = Y1;
	}

	// adjust the order

	if (X2 < X1)
	{
		double c = X1;
		X1 = X2;
		X2 = c;
	}
	if (Y2 < Y1)
	{
		double c = Y1;
		Y1 = Y2;
		Y2 = c;
	}

	//// enter the ROI
	double weight = 0, Ray = 0;
	short GridX = 0, GridY = 0;        // candidate crosspoint index in matlab(0~256)
	double GridY_x = 0, GridX_y = 0;    // candidate crosspoint index in matlab(0~256)
	short DetectPoint_x = 0, DetectPoint_y = 0; // current pixel index in matlab pixel index in matlab(0~255)
	int	Pointid = 0;   
	double XCross = XStart / dev_resolution[1], YCross = YStart / dev_resolution[0];     // current crosspoint index in matlab(0~256)

	//while (((XCross * dev_resolution[1]) >= Range[0]) && ((XCross * dev_resolution[1]) <= Range[2]) 
	//	&& ((YCross * dev_resolution[0]) >= Range[1]) && ((YCross * dev_resolution[0]) <= Range[3]))
	for (short i = 0; i < (height + width - 1); i++)
	{
		//dev_Projection[threadid] = 10000;

		// judge whether XCross/YCross is integer
		if (XCross - (double)((short)XCross) < EPS)
		{
			GridX = XCross + x_signal;			
		}
		else
		{
			GridX = floor(XCross) + flag1to1or_1to0(x_signal);
		}
		GridX_y = (DetectPoint_ystart + (GridX * dev_resolution[1] - DetectPoint_xstart) * X2Y  /** tan(theta + PI / 2)*/) / dev_resolution[0];

		if (YCross - (double)((short)YCross) < EPS)
		{
			GridY = YCross + y_signal;			
		}
		else
		{
			GridY = floor(YCross) + flag1to1or_1to0(y_signal);
		}
		GridY_x = (DetectPoint_xstart + (GridY * dev_resolution[0] - DetectPoint_ystart) * Y2X /*/ (tan(theta + PI / 2) + EPS)*/) / dev_resolution[1];

		//judge which crosspoint is the nearest
		if (Distancesq(GridX, GridX_y, XCross, YCross) >= Distancesq(GridY_x, GridY, XCross, YCross))
		{
			weight = sqrt(Distancesq(GridY_x * dev_resolution[1], GridY * dev_resolution[0], XCross * dev_resolution[1], YCross * dev_resolution[0]));
			DetectPoint_x = floor(MID(GridY_x, XCross));                 // the midpoint locates the pixel
			DetectPoint_y = floor(MID(GridY, YCross));
			XCross = GridY_x;    // update
			YCross = GridY;
		}
		else
		{
			weight = sqrt(Distancesq(GridX * dev_resolution[1], GridX_y * dev_resolution[0], XCross * dev_resolution[1], YCross * dev_resolution[0]));
			DetectPoint_x = floor(MID(GridX, XCross));           
			DetectPoint_y = floor(MID(GridX_y, YCross));
			XCross = GridX;   // update
			YCross = GridX_y;
		}

		//judge whether the point is in the ROI
		if ((DetectPoint_x >= 0) && (DetectPoint_x < width) && (DetectPoint_y >= 0) && (DetectPoint_y < height ))
		{
			Pointid = DetectPoint_x * height + DetectPoint_y;
			Ray += weight * dev_Pic[Pointid];
			// for debug{
			//dev_Projection[threadid] = 8000;
			//break;
			// }
		}
		else
		{
			// for debug{
			//dev_Projection[threadid] = 9000;
			//}
			break;
		}
			
	}

	//__syncthreads();
	//for debug{
	//dev_Projection[threadid] = Ray;
	//}

	dev_Projection[threadid] = Ray;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t ProjectionParallel_2D(const float *Pic, double *Projection, const float *thetaRange, const float *t_Range,
	const int height, const int width, const double Center_y, const double Center_x, const int Ltheta, const int Lt,
	const double *resolution)
{
	mexPrintf("Hello GenMatParalell!\n");
	float *dev_Pic = 0, * dev_thetaRange = 0, *dev_t_Range = 0 ;
	double *dev_Projection = 0, *dev_resolution = 0;

	int threadcubic_x = MIN(threadX, Lt);
	int blockcubic_x = MIN(blockX, Ltheta);
	int LtResidual = Lt % threadX;
	int LthetaResidual = Ltheta % blockX;
	int tstart = 0;
	int thetastart = 0;
	int tTime = Lt / threadX;
	int thetaTime = Ltheta / blockX;

	const dim3 thread_cubic(threadcubic_x, 1, 1);
	const dim3 block_cubic(blockcubic_x, 1, 1);
	dim3 thread_cubic_residual(1, 1, 1);  // initial
	dim3 block_cubic_residual(1, 1, 1);  // initial

	mexPrintf("%d %d %d %d\n", threadcubic_x, blockcubic_x, LtResidual, LthetaResidual);

	if (LtResidual != 0)
	{
		thread_cubic_residual.x = LtResidual;
	}
	if (LthetaResidual != 0)
	{
		block_cubic_residual.x = LthetaResidual;
	}

	cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		mexPrintf("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
        goto Error;
    }

	mexPrintf("Call for GPU space.\n");

    // Allocate GPU buffers for three vectors (two input, one output).
	cudaStatus = cudaMalloc((void**)&dev_Pic, height * width * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_Pic cudaMalloc failed!");
		mexPrintf("dev_Pic cudaMalloc failed!\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_Projection, Ltheta * Lt * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_Projection cudaMalloc failed!");
		mexPrintf("dev_Projection cudaMalloc failed!\n");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_thetaRange, Ltheta * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_thetaRange cudaMalloc failed!");
		mexPrintf("dev_thetaRange cudaMalloc failed!\n");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_t_Range, Lt * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "dev_t_Range cudaMalloc failed!");
		mexPrintf("dev_t_Range cudaMalloc failed!\n");
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
	cudaStatus = cudaMemcpy(dev_Pic, Pic, height * width * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "thetaRange cudaMemcpy failed!");
		mexPrintf("thetaRange cudaMemcpy failed!\n");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_thetaRange, thetaRange, Ltheta * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "thetaRange cudaMemcpy failed!");
		mexPrintf("thetaRange cudaMemcpy failed!\n");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_t_Range, t_Range, Lt * sizeof(float), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "t_Range cudaMemcpy failed!");
		mexPrintf("t_Range cudaMemcpy failed!\n");
        goto Error;
    }

	cudaStatus = cudaMemcpy(dev_resolution, resolution, 2 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "dev_resolution cudaMemcpy failed!");
		mexPrintf("dev_resolution cudaMemcpy failed!\n");
		goto Error;
	}

	mexPrintf("Launch computation projection of each lines.\n");
    // Launch a kernel on the GPU with one thread for each element.
	
	for (int numt = 0; numt < tTime; numt++)
	{
		for (int numtheta = 0; numtheta < thetaTime; numtheta++)
		{
			tstart = numt * threadX;
			thetastart = numtheta * blockX;
			//mexPrintf("%d %d\n", Pstart, Xigamastart);
			ProjectionParallel << <block_cubic, thread_cubic >> > (dev_Pic, dev_Projection, dev_t_Range, dev_thetaRange,
				Center_y, Center_x, dev_resolution, height, width, Lt, tstart, thetastart);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionCone launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionCone launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at tstarttstart: %d thetastart: %d \n",
					tstart, thetastart);
				goto Error;
			}
		}
	}

	if (LtResidual != 0)
	{
		tstart = Lt - LtResidual;
		if (LthetaResidual != 0)
		{
			thetastart = Ltheta - LthetaResidual;
			//("%d %d\n", Pstart, Xigamastart);
			ProjectionParallel << <block_cubic_residual, thread_cubic_residual >> > (dev_Pic, dev_Projection, dev_t_Range, dev_thetaRange,
				Center_y, Center_x, dev_resolution, height, width, Lt, tstart, thetastart);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionCone launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionCone launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at tstarttstart: %d thetastart: %d \n",
					tstart, thetastart);
				goto Error;
			}
		}

		for (int numtheta = 0; numtheta < thetaTime; numtheta++)
		{
			thetastart = numtheta * blockX;
			//("%d %d\n", Pstart, Xigamastart);
			ProjectionParallel << <block_cubic, thread_cubic_residual >> > (dev_Pic, dev_Projection, dev_t_Range, dev_thetaRange,
				Center_y, Center_x, dev_resolution, height, width, Lt, tstart, thetastart);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionCone launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionCone launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at tstarttstart: %d thetastart: %d \n",
					tstart, thetastart);
				goto Error;
			}
		}
	}
	if (LthetaResidual != 0)
	{
		thetastart = Ltheta - LthetaResidual;
		for (int numt = 0; numt < tTime; numt++)
		{
			tstart = numt * threadX;
			//mexPrintf("%d %d\n", Pstart, Xigamastart);
			ProjectionParallel << <block_cubic_residual, thread_cubic >> > (dev_Pic, dev_Projection, dev_t_Range, dev_thetaRange,
				Center_y, Center_x, dev_resolution, height, width, Lt, tstart, thetastart);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "ProjectionCone launch failed: %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("ProjectionCone launch failed %s\n", cudaGetErrorString(cudaStatus));
				mexPrintf("Error happens at tstarttstart: %d thetastart: %d \n",
					tstart, thetastart);
				goto Error;
			}
		}
	}

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Projection launch failed: %s\n", cudaGetErrorString(cudaStatus));
		mexPrintf("Projection launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
   
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		mexPrintf("cudaDeviceSynchronize returned error\n");
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(Projection, dev_Projection, Ltheta * Lt * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!\n");
        goto Error;
    }

Error:
    cudaFree(dev_t_Range);
    cudaFree(dev_thetaRange);
    cudaFree(dev_Projection);
	cudaFree(dev_Pic);
	cudaFree(dev_resolution);

    return cudaStatus;
}
