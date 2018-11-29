#pragma once
#pragma once
// 2018/11/16 apply GPU acceleration
#ifndef __BACKPROJECTION_H__
#define __BACKPROJECTION_H__
#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#include "G:\CUDA\Development\include\cuda_runtime.h"
//#include "G:\CUDA\Development\include\device_launch_parameters.h"
#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\include\cuda_runtime.h"
#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\include\device_launch_parameters.h"

//  Maximum threads of each dimension of a block: 1024 x 1024 x 64
// Maximum threads of each dimension of a grid: 2147483647 x 65535 x 65535
// Maximum threads of each dimension of a grid: 2097152(1024*2048) x 64 x 1024
using namespace std;

#define threadX 512
#define blockX 512

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y
#define pow2(x) (1.0*(x)*(x))
#define Distance(x1,y1,z1,x2,y2,z2) (sqrt(pow2(x1-x2)+pow2(y1-y2)+pow2(z1-z2)))

cudaError_t BackPro(float *Display, const float *R, const float *Xigamadomain, const float *Pdomain,
	const float *BetaScanRange, const double Distance, const int LBeta, const int LP, const int LXigama,
	const double *Size, const int t_length, const int s_length, const int z_length);
#endif
