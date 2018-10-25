#pragma once
// 2018/04/16
#ifndef __PROJECTIONCONE_3D_H__
#define __PROJECTIONCONE_3D_H__
#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "G:\CUDA\Development\include\cuda_runtime.h"
#include "G:\CUDA\Development\include\device_launch_parameters.h"
using namespace std;
#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)
#define flag1to1or_1to0(flag) ((1 + flag)/2)    //when flag=1, ans=1; when flag=-1, ans= 0
#define pow2(x) (1.0*(x)*(x))
#define ABS(x) (x>0?x:(-x))
#define Distancesq(x1,y1,z1,x2,y2,z2) (pow2(x1-x2)+pow2(y1-y2)+pow2(z1-z2))
#define MID(x,y) ((x+y)/2)
#define threadX 512
#define blockX 512

cudaError_t ProjectionCone_3D(const float *Pic, float *Projection, const float *BetaScanRange, const float *Pdomain,
	const float *Xigamadomain, const int t_length, const int s_length, const int z_length, const double Center_t,
	const double Center_s, const double Center_z, const int LBeta, const int LP, const int LXigama, const double Distance,
	const double *resolution);

#endif