#pragma once
#ifndef __IMPLEMENTUI_H__
#define __IMPLEMENTUI_H__
#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\include\cuda_runtime.h"
#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\include\device_launch_parameters.h"
//  Maximum threads of each dimension of a block: 1024 x 1024 x 64
// Maximum threads of each dimension of a grid: 2147483647 x 65535 x 65535
// Maximum threads of each dimension of a grid: 2097152(1024*2048) x 64 x 1024
using namespace std;
#define THREAD_SIZE_X 936
#define Filter_SIZE (2*THREAD_SIZE_X-1)
//#define THREAD_SIZE_Y 4
//#define THREAD_SIZE_Z 32
#define BLOCK_SIZE_X 80
#define BLOCK_SIZE_Y 800
//#define BLOCK_SIZE_Z 25
#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y
//#define pow2(x) pow(x,2)

//#define ARRAY_SIZE_IN_BYTES ((THREAD_SIZE_X*BLOCK_SIZE_X)*(THREAD_SIZE_Y*BLOCK_SIZE_Y)*(THREAD_SIZE_Z*BLOCK_SIZE_Z)*sizeof (double)) 
//#define ARRAY_SIZE_OUT_BYTES ((THREAD_SIZE_X*BLOCK_SIZE_X)*(THREAD_SIZE_Y*BLOCK_SIZE_Y)*(THREAD_SIZE_Z*BLOCK_SIZE_Z)*sizeof (double)) 
#define ARRAY_SIZE_IN_BYTES ( THREAD_SIZE_X * BLOCK_SIZE_X * BLOCK_SIZE_Y * sizeof (double)) 
#define ARRAY_SIZE_OUT_BYTES ( THREAD_SIZE_X * BLOCK_SIZE_X * BLOCK_SIZE_Y * sizeof (double)) 

cudaError_t FDKUI(double *Display, const double *R, const double *Xigamadomain, const double *Gamadomain,
	const double *BetaScanRange, const double betaStartAngle, const double MaxGama, const double Distance,
	const double Distance_s2d, const int LGama, const int LBeta, const int LXigama, const int *Size, const int z_length);

#endif