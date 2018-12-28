#include "Backprojection.h"
// 2018/11/16 apply GPU accelerationd
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

	//mexPrintf("hello!\n");
	float *R = (float*)mxGetPr(prhs[0]);
	const float *Xigamadomain = (float*)mxGetPr(prhs[1]);
	const int LXigama = (int)mxGetM(prhs[1]);
	const float *Pdomain = (float*)mxGetPr(prhs[2]);
	const int LP = (int)mxGetM(prhs[2]);
	const float *BetaScanRange = (float*)mxGetPr(prhs[3]);
	const int LBeta = (int)mxGetM(prhs[3]);
	const double Distance = mxGetScalar(prhs[4]);
	const double *Size = mxGetPr(prhs[5]);      // Size denotes the size of Display matrix(3d)
	const int t_length = (int)mxGetScalar(prhs[6]);
	const int s_length = (int)mxGetScalar(prhs[7]);
	const int z_length = (int)mxGetScalar(prhs[8]);

	//mexPrintf("%f \n", Distance);
	//mexPrintf("%d %d %d\n", t_length, s_length ,z_length);
	/*mexPrintf("%f\n", R[9781670]);*/

	const long LDisplay = t_length * s_length * z_length;
	const long LR = LBeta * LP * LXigama;

	float *Displaycuda = new float[LDisplay];
	memset(Displaycuda, 0, sizeof(float));

	plhs[0] = mxCreateDoubleMatrix(LDisplay, 1, mxREAL);
	double *Display = mxGetPr(plhs[0]);

	cudaError_t cudaStatus = BackPro(Displaycuda, R, Xigamadomain, Pdomain, BetaScanRange, Distance, LBeta,
		LP, LXigama, Size, t_length, s_length, z_length);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BackPro failed!\n");
		mexPrintf("BackPro failed!\n");
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
}