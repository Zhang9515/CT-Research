#pragma once
// 2018/04/09
#ifndef __PROJECTIONPARALLEL_2D_H__
#define __PROJECTIONPARALLEL_2D_H__
#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
// if local environment variable has been set, there is no need to give the complete including path 
// lab computer
//#include "G:\CUDA\Development\include\cuda_runtime.h"
//#include "G:\CUDA\Development\include\device_launch_parameters.h"
// server2(maybe has changed)
/*C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\include\*/
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;

#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)
#define flag1to1or_1to0(flag) ((1 + flag)/2)    //when flag=1, ans=1; when flag=-1, ans= 0
//#define flag1to1or_1to0(flag)((-1 + flag)/2)    //when flag=1, ans=0; when flag=-1, ans= -1
#define pow2(x) (1.0*(x)*(x))
#define ABS(x) (x>0?x:(-x))
#define Distancesq(x1,y1,x2,y2) (pow2(x1-x2)+pow2(y1-y2))
#define MID(x,y) ((x+y)/2)
#define threadX 512
#define blockX 512
#define Maxlim(x) (x>1e+16?1e+16:x)
#define Minlim(x) (x<-1e+16?-1e+16:x)

cudaError_t ProjectionParallel_2D(const float *Pic, double *Projection, const double *thetaRange, const double *t_Range,
	const int height, const int width, const double Center_y, const double Center_x, const int Ltheta, const int Lt,
	const double *resolution);

#endif