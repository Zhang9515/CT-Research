#pragma once
#pragma once
// 2019/03/28 apply GPU acceleration
#ifndef __FBPPARALLEL_H__
#define __FBPPARALLEL_H__
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
//  Maximum threads of each dimension of a block: 1024 x 1024 x 64
// Maximum threads of each dimension of a grid: 2147483647 x 65535 x 65535
// Maximum threads of each dimension of a grid: 2097152(1024*2048) x 64 x 1024
using namespace std;

#define threadX 256
#define blockX 256
#define Filterlengthlimit 2048

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y

cudaError_t FBPpara(float *Display, const float *R, const float *Tdomain, const float *ThetaScanRange, const int LTheta,
	const int LT, const double *Size, const int t_length, const int s_length);

#endif