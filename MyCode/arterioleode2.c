#include "mex.h"
#include "math.h"


void arterioleODE2()
{






	outPut
	
  
}



// nlhs -> Number of output variables
// plhs -> Array of mxArray pointers to output variables
// nrhs -> Number of input variables
// prhs -> Array of mxArray pointers to the input variables


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {


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


plhs[0] = mxCreateDoubleMatrix(2,1, mxREAL);
outPut = mxGetPr(plhs[0]);

arterioleODE2(t,u,Pin,Qnum,Pdrop,multiplier,capdensity,CMD_DA,td,Pac,... 
Pend,MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,Pdistc,Params_dist,CLc,ta, outPut);

}
