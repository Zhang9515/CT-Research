#pragma once
// 2018.8.6
#ifndef __ITERATIVEFDKERRCOMPUTE_H__
#define __ITERATIVEFDKERRCOMPUTE_H__

#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "G:\CUDA\Development\include\cuda_runtime.h"
#include "G:\CUDA\Development\include\device_launch_parameters.h"

#define Round(x) (x>0)?floor(x+0.5):ceil(x-0.5)
#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)
#define pow2(x) (1.0*(x)*(x))
#define threadX 512
#define blockX 512
#define ABS(x) ((x)>0?(x):(-(x)))

using namespace std;

cudaError_t IterativeFDKerrCompute(const float *Pic, float *ErrorSlicecuda, const float *BetaScanRange, const float *Udomain,
	const double *Pdomain, const double *G, const double *Size, const int t_length, const int s_length, const int z_length,
	const double Center_t, const double Center_s, const double Center_z, const int LBeta, const int LP, const int LU,
	const int LG, const double Distance, const double *resolution2, const int t, const float BetaScanInt, const float dU, 
	const double PInt);

#endif