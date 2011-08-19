#include "mex.h"
#include "config.h"
#include "ltfat.h"

/* Calling convention:
 *  phase=comp_heapint(s,itime,ifreq,a,tol);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int a, M, N, L, W;
   double tol;
   
   double *s, *tgrad, *fgrad, *phase;
   
      /* Get inputs */
   s     = mxGetPr(prhs[0]);
   tgrad = mxGetPr(prhs[1]);
   fgrad = mxGetPr(prhs[2]);
   a     = (int)mxGetScalar(prhs[3]);
   tol   = mxGetScalar(prhs[4]);

   /* Get matrix dimensions.*/
   M = mxGetM(prhs[0]); 
   N = mxGetN(prhs[0]);
   L = N*a;
   W = 1;

   /* Create output matrix and zero it.*/
   plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);   

   /* Get pointer to output */
   phase=mxGetPr(plhs[0]);
   
   heapint((const double*)s,
	   (const double*)tgrad,
	   (const double*)fgrad,
	   a, M, L, W, 
	   tol, phase);

   return;
   
}


