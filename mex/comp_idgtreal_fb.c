#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"


/* Calling convention:
 *  comp_idgtreal_fb(coef,g,L,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, W, a, M, N, M2, gl;
   ltfat_complex *c_combined;
   
   int ii;
   
   /* Get matrix dimensions.*/
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetN(prhs[0])/N;
   gl = mxGetM(prhs[1]); 


   M2= M/2+1;
   
   /* Create temporary matrices to convert to correct complex layout. */
   c_combined=mxMalloc(M2*N*W*sizeof(ltfat_complex));   

   split2combined(M2*N*W, prhs[0], c_combined);

   plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);
   
   idgtreal_fb((const ltfat_complex*)c_combined,
		(const double*)mxGetPr(prhs[1]),L,gl,W,a,M,
		(double*)mxGetPr(plhs[0]));
   
   mxFree(c_combined);
   
   return;
   
}


