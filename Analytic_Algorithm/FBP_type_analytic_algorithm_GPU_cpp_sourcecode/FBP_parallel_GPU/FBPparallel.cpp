#include "FBPparallel.h"
// 2019/03/28 apply GPU accelerationd
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

	//mexPrintf("hello!\n");
	float *R = (float*)mxGetPr(prhs[0]);
	const float *ThetaScanRange = (float*)mxGetPr(prhs[1]);
	const int LTheta = (int)mxGetM(prhs[1]);
	const float *Tdomain = (float*)mxGetPr(prhs[2]);
	const int LT = (int)mxGetM(prhs[2]);
	const double *Size = mxGetPr(prhs[3]);      // Size denotes the actual size of Display matrix(2d)
	const int t_length = (int)mxGetScalar(prhs[4]);
	const int s_length = (int)mxGetScalar(prhs[5]);

	//mexPrintf("size:%d\n", (int)mxGetM(prhs[3]));
	//mexPrintf("%f \n", Distance);
	//mexPrintf("%d %d %d\n", t_length, s_length ,z_length);
	/*mexPrintf("%f\n", R[9781670]);*/

	const long LDisplay = t_length * s_length;
	const long LR = LTheta * LT;

	float *Displaycuda = new float[LDisplay];
	memset(Displaycuda, 0, sizeof(float));

	plhs[0] = mxCreateDoubleMatrix(LDisplay, 1, mxREAL);
	double *Display = mxGetPr(plhs[0]);
	//for debug{
	/*float *Rcovcuda = new float[LR];
	memset(Rcovcuda, 0, sizeof(float));*/

	//plhs[0] = mxCreateDoubleMatrix(LR, 1, mxREAL);
	//double *Rcov = mxGetPr(plhs[0]);
	//}end debug

	cudaError_t cudaStatus = FBPpara(Displaycuda, R, Tdomain, ThetaScanRange, LTheta,
		LT, Size, t_length, s_length);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FBP failed!\n");
		mexPrintf("FBP failed!\n");
		return;
	}

	// mxCreateDoubleMatrix: store unit is double, every output should be converted to double explicitly 

	for (long num = 0; num < LDisplay; num++)
		Display[num] = (double)Displaycuda[num];

	//memcpy(Display, Displaycuda, LDisplay * sizeof(float));

	//mexPrintf("%lf\n", RPreweight[0]);

	//cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!\n");
		mexPrintf("cudaDeviceReset failed!\n");
		return;
	}
	delete[] Displaycuda;
	/*delete[] Rcovcuda;*/
}