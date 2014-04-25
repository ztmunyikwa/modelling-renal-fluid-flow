#include "mex.h"
#include "math.h"


void arterioleODE2(double y, double t, double *z)
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

int main(int argc, char* argv[]) {
//argc - number of command line arguments
//argv - the command line arguments as an array 

//function z = arterioleODE(t,u,Pin,Qnum,Pdrop,multiplier,capdensity,CMD_DA,td,Pac,... 
	//Pend,MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,Pdistc,Params_dist,CLc,ta)

if (MyoOn_dist == 1) {
	Pstartdist = Pdistnew;
} else {

	Pstartdist = Pdistc;
	exit(EXIT_SUCCESS);
}

if (TGFOn_dist == 1) {
	CMD_DA =CMD_DA;

} else {
	CMD_DA = CLc;
	exit(EXIT_SUCCESS);
}



}

//MATLAB Code
// if MyoOn_dist == 1
//     Pstartdist = Pdistnew;
// else
//     Pstartdist = Pdistc;  %control state assumption
// %     Pstartdist = Pdistnew;  %completely zeroed out
// %     Params_dist(6) = 0;     %completely zeroed out
// end

// if TGFOn_dist == 1
//     CMD_DA = CMD_DA;
// else
//     CMD_DA = CLc;   %control state assumption
// %     Params_dist(7) = 0;  %completely zero out
// end
