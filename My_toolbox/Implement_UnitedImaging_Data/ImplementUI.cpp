#include "ImplementUI.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

	//mexPrintf("hello!\n");
	double *R = mxGetPr(prhs[0]);
	const int RSize = mxGetM(prhs[0]);
	const double *Xigamadomain = mxGetPr(prhs[1]);
	const int LXigama = mxGetM(prhs[1]);
	const double *Gamadomain = mxGetPr(prhs[2]);
	const int LGama = mxGetM(prhs[2]);
	const double *BetaScanRange = mxGetPr(prhs[3]);
	const int LBeta = mxGetM(prhs[3]);
	const double betaStartAngle = mxGetScalar(prhs[4]);
	const double MaxGama = mxGetScalar(prhs[5]);
	const double Distance = mxGetScalar(prhs[6]);
	const double Distance_s2d = mxGetScalar(prhs[7]);
	const double *Size_d = mxGetPr(prhs[8]);      // Size denotes the size of Display matrix(3d)
	int *Size = new int[3];
	for (int i = 0; i < 3; i++)
		Size[i] = (int) Size_d[i];
	const int z_length = (int) mxGetScalar(prhs[9]);

	const long LDisplay = Size[0] * Size[1] * z_length;

	/*mexPrintf("%d %lf %lf %lf %lf %lf %lf\n", RSize , Xigamadomain[1] , Gamadomain[1] , BetaScanRange[1] , betaStartAngle , 
		MaxGama , Distance);*/

	double *Displaycuda = new double[LDisplay];
	memset(Displaycuda, 0, sizeof(double));

	plhs[0] = mxCreateDoubleMatrix(LDisplay, 1, mxREAL);
	double *Display = mxGetPr(plhs[0]);

	cudaError_t cudaStatus = FDKUI(Displaycuda, R, Xigamadomain, Gamadomain,
		BetaScanRange, betaStartAngle, MaxGama, Distance, Distance_s2d, LGama, LBeta, LXigama, Size, z_length);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FDKUI failed!\n");
		mexPrintf("FDKUI failed!\n");
	    return;
	}

	memcpy(Display, Displaycuda, LDisplay * sizeof(double));

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