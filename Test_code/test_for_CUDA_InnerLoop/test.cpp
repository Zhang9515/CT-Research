#include "test.h"
// 2018/04/16
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
		const float *R = (float*)mxGetPr(prhs[0]);
		
		float *Displaycuda = new float[256];
		memset(Displaycuda, 0, sizeof(float));

		//mexPrintf("%d %d %d\n", t_length, s_length, s_length);
		//mexPrintf("%f %f %f %d %f %d %f %d %ld %f\n", Pic[77617], Size[2], BetaScanRange[0], LBeta, Pdomain[0], 
		//	LP, Xigamadomain[0], LXigama, LMat, Distance);

		plhs[0] = mxCreateDoubleMatrix(256, 1, mxREAL);
		double *Display = mxGetPr(plhs[0]);

		// Add vectors in parallel.
		cudaError_t cudaStatus = FDKpro(Displaycuda, R);

		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "GenMatParalell failed!");
			mexPrintf("GenMatParalell failed!\n");
			return;
		}

		// mxCreateDoubleMatrix: store unit is double, every output should be converted to double explicitly 

		for (int num = 0; num < 256; num++)
			Display[num] = (double)Displaycuda[num];

		//memcpy(ProjectionCone, ProjectionConecuda, LMat * sizeof(double));

		// cudaDeviceReset must be called before exiting in order for profiling and
		// tracing tools such as Nsight and Visual Profiler to show complete traces.
		cudaStatus = cudaDeviceReset();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceReset failed!\n");
			return;
		}
		delete [] Displaycuda;
}