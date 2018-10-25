#include "ProjectionParallel_2D.h"
// 2018/04/09
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
		const double *Pic = mxGetPr(prhs[0]);
		const int height = (int)mxGetScalar(prhs[1]);
		const int width = (int)mxGetScalar(prhs[2]);
		const double *Size = mxGetPr(prhs[3]);
		const double Center_y = Size[0]/2;    // vertical
		const double Center_x = Size[1]/2;    // horizontal
		const double *thetaRange = mxGetPr(prhs[4]);
		const int Ltheta = mxGetM(prhs[4]);
		const double *t_range = mxGetPr(prhs[5]);
		const int Lt = mxGetM(prhs[5]);
		const long LMat = Ltheta * Lt;

		const double resolution[] = {Size[0] / height, Size[1] / width};
		double *ProjectionParallelcuda = new double[LMat];
		memset(ProjectionParallelcuda, 0, sizeof(double));

		//mexPrintf("%lf %lf %ld\n", Center_y,Center_x, LMat);

		plhs[0] = mxCreateDoubleMatrix(LMat, 1, mxREAL);
		double *ProjectionParallel = mxGetPr(plhs[0]);

		// Add vectors in parallel.
		cudaError_t cudaStatus = ProjectionParallel_2D(Pic, ProjectionParallelcuda, thetaRange,
			t_range, height, width, Center_y, Center_x, Ltheta, Lt, resolution);

		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "GenMatParalell failed!");
			mexPrintf("GenMatParalell failed!\n");
			return;
		}

		memcpy(ProjectionParallel, ProjectionParallelcuda, LMat * sizeof(double));

		// cudaDeviceReset must be called before exiting in order for profiling and
		// tracing tools such as Nsight and Visual Profiler to show complete traces.
		cudaStatus = cudaDeviceReset();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceReset failed!\n");
			return;
		}
		delete [] ProjectionParallelcuda;
}