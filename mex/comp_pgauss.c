#include "mex.h"
#include "config.h"
#include "ltfat.h"

/* Calling convention:
 *  comp_pgauss(L,w,center);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int L;
   double w, center;
   double *g;
  
   L=(int)mxGetScalar(prhs[0]);
   w=(double)mxGetScalar(prhs[1]);
   center=(double)mxGetScalar(prhs[2]);

   plhs[0] = mxCreateDoubleMatrix(L, 1, mxREAL);
   
   g = mxGetPr(plhs[0]);
  
   pgauss(L, w, center,(double*)g);

   return;
  
}


