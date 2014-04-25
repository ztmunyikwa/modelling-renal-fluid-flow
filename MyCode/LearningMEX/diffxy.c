//To get a solution, we can type:
// [T, XY] ode45('diffxy', 0, 10, [0 1 0])


void diffxy(double x, double *y, double *z, mwSize n)
{
  mwSize i;
  
  for (i=0; i<n; i++) {
    z[i] = x * y[i];
  }
}



// nlhs -> Number of output variables
// plhs -> Array of mxArray pointers to output variables
// nrhs -> Number of input variables
// prhs -> Array of mxArray pointers to the input variables


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

}