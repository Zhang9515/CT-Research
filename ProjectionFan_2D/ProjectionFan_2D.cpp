#include "ProjectionFan_2D.h"
// 2018/04/16
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	const float *Pic = (float*)mxGetPr(prhs[0]);
	const int t_length = (int)mxGetScalar(prhs[1]);
	const int s_length = (int)mxGetScalar(prhs[2]);
	const double *Size = mxGetPr(prhs[4]);
	const double Center_t = Size[0] / 2;    // horizontal plane
	const double Center_s = Size[1] / 2;    // horizontal plane
	const double *BetaScanRange = (double*)mxGetPr(prhs[5]);
	const int LBeta = mxGetM(prhs[5]);
	const double *Pdomain = (double*)mxGetPr(prhs[6]);
	const int LP = mxGetM(prhs[6]);
	const long LMat = LBeta * LP;
	const double Distance = (double)mxGetScalar(prhs[7]);

	const double resolution[] = { Size[0] / t_length, Size[1] / s_length };
	double *ProjectionConecuda = new double[LMat];
	memset(ProjectionConecuda, 0, sizeof(double));

	//mexPrintf("%d %d %d\n", t_length, s_length, s_length);
	//mexPrintf("%f %f %f %d %f %d %f %d %ld %f\n", Pic[77617], Size[2], BetaScanRange[0], LBeta, Pdomain[0], 
	//	LP, Xigamadomain[0], LXigama, LMat, Distance);

	plhs[0] = mxCreateDoubleMatrix(LMat, 1, mxREAL);
	double *ProjectionCone = mxGetPr(plhs[0]);

	// Add vectors in parallel.
	cudaError_t cudaStatus = ProjectionCone_3D(Pic, ProjectionConecuda, BetaScanRange,
		Pdomain, t_length, s_length, Center_t, Center_s, LBeta, LP, Distance, resolution);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "GenMatParalell failed!");
		mexPrintf("GenMatParalell failed!\n");
		return;
	}

	// mxCreateDoubleMatrix: store unit is double, every output should be converted to double explicitly 

	for (long num = 0; num < LMat; num++)
		ProjectionCone[num] = (double)ProjectionConecuda[num];

	//memcpy(ProjectionCone, ProjectionConecuda, LMat * sizeof(double));

	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tooLP such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!\n");
		return;
	}
	delete[] ProjectionConecuda;
}