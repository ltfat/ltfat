#include "mex.h"
#include "config.h"
#include "ltfat.h"

/* Calling convention:
 *  cout=comp_gabreassign(s,itime,ifreq,a);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int a, M, N, L;
   double *s,*tgrad, *fgrad,*sr;
   
   int ii;
   
   /* Get matrix dimensions.*/
   M = mxGetM(prhs[0]); 
   N = mxGetN(prhs[0]);
   a = (int)mxGetScalar(prhs[3]);
   L = N*a;

   s     = mxGetPr(prhs[0]);   
   tgrad = mxGetPr(prhs[1]);   
   fgrad = mxGetPr(prhs[2]);   
   
   plhs[0] = mxCreateDoubleMatrix(M,N, mxREAL);
   sr      = mxGetPr(plhs[0]);   

   gabreassign((double*)s,(double*)tgrad,(double*)fgrad,
               L,1,a,M,(double*)sr);
                  
   return;
  
}


