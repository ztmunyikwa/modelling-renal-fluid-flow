
#include "mex.h"
#include "math.h"

void fun1(double y, double t, double *z)
{
	double temp = 2 - pow(y, 2);
	
	z[0] = -t*y/sqrt(temp); 

  
}



// nlhs -> Number of output variables
// plhs -> Array of mxArray pointers to output variables
// nrhs -> Number of input variables
// prhs -> Array of mxArray pointers to the input variables


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
	

	double y;
	double t;
	double z;
	double * outPut;

	if(nrhs !=2) {
		mexErrMsgIdAndTxt("MyToolBox:fun1:nrhs", "Two inputs required");
	}
	if(nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolBox:fun1:nlhs","One output required");
	}

	/*make sure inputs are both scalars */

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("MyToolbox:fun1:notScalar",
                      "Inputs must be a scalars.");
	}

	t = mxGetScalar(prhs[0]);

	y = mxGetScalar(prhs[1]);

	/*create the output array */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

	outPut = mxGetPr(plhs[0]);

	fun1(y, t, outPut);


}