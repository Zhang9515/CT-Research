#pragma once
// 2018/04/16
#ifndef __PROJECTIONFAN_2D_H__
#define __PROJECTIONFAN_2D_H__
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
#define pow2(x) (1.0*(x)*(x))
#define ABS(x) (x>0?x:(-x))
#define Distancesq(x1,y1,x2,y2) (pow2(x1-x2)+pow2(y1-y2))
#define MID(x,y) ((x+y)/2)
#define threadX 256
#define blockX 256
#define Filterlengthlimit 900
#define Maxlim(x) (x>1e+16?1e+16:x)
#define Minlim(x) (x<-1e+16?-1e+16:x)

cudaError_t ProjectionCone_3D(const float *Pic, double *Projection, const double *BetaScanRange, const double *Pdomain,
	const int t_length, const int s_length, const double Center_t, const double Center_s, const int LBeta, const int LP, 
	const double Distance, const double *resolution);

#endif