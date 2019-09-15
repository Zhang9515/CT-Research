#include "IterativeFDKerrCompute.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
		const float *Pic = (float*)mxGetPr(prhs[0]);
		const int t_length = (int)mxGetScalar(prhs[1]);
		const int s_length = (int)mxGetScalar(prhs[2]);
		const int z_length = (int)mxGetScalar(prhs[3]);
		const double *Size = mxGetPr(prhs[4]);
		const double Center_t = Size[0] / 2;    // horizontal plane
		const double Center_s = Size[1] / 2;    // horizontal plane
		const double Center_z = Size[2] / 2;    // horizontal plane
		const float *BetaScanRange = (float*)mxGetPr(prhs[5]);
		const int LBeta = mxGetM(prhs[5]);
		const float BetaScanInt = ABS(BetaScanRange[1] - BetaScanRange[0]);
		const float *Udomain = (float*)mxGetPr(prhs[6]);
		const int LU = mxGetM(prhs[6]);
		const float dU = ABS(Udomain[1] - Udomain[0]);
		const double *Pdomain = mxGetPr(prhs[7]);
		const int LP = mxGetM(prhs[7]);
		const double PInt = ABS(Pdomain[1] - Pdomain[0]);
		const long LMat = LBeta * s_length * z_length;
		const double Distance = mxGetScalar(prhs[8]);
		const double *G = mxGetPr(prhs[9]);
		const int LG = mxGetM(prhs[9]);
		const double Resolution2[] = {Size[0] / t_length, Size[1] / s_length, Size[2] / z_length };
		const int t = mxGetScalar(prhs[10]) - 1;    // input index is matlab index; start from 1 different from C++ 0

		float *ErrorSlicecuda = new float[LMat];
		memset(ErrorSlicecuda, 0, sizeof(float));

		plhs[0] = mxCreateDoubleMatrix(LMat, 1, mxREAL);
		double *ErrorSlice = mxGetPr(plhs[0]);
		
		 // GPU-CPU interface
		cudaError_t cudaStatus = IterativeFDKerrCompute(Pic, ErrorSlicecuda, BetaScanRange,
			Udomain, Pdomain, G, Size, t_length, s_length, z_length, Center_t, Center_s, Center_z,
			LBeta, LP, LU, LG,  Distance, Resolution2, t, BetaScanInt, dU, PInt);

		mexPrintf("Parallel process end!\n");

		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "ParalellCompute failed!");
			mexPrintf("ParalellCompute failed!\n");
			return;
		}

		// convert the data format for matlab
		for (long num = 0; num < LMat; num++)
			ErrorSlice[num] = (double)ErrorSlicecuda[num];

		cudaStatus = cudaDeviceReset();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceReset failed!");
			return;
		}
}