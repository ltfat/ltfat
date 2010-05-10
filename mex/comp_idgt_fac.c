#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_idgt_fac(coef,gf,L,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, W, a, M, N;
   ltfat_complex *f_combined, *gf_combined, *c_combined;
   
   /* Get matrix dimensions.*/
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetM(prhs[0])*mxGetN(prhs[0])/(M*N); 
   
   /* Create temporary matrices to convert to correct complex layout. */

   c_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));   
   gf_combined=mxMalloc(L*sizeof(ltfat_complex));
   f_combined=mxMalloc(L*W*sizeof(ltfat_complex));
   
   split2combined(M*N*W, prhs[0], c_combined);
   split2combined(L, prhs[1], gf_combined);
   
   idgt_fac((const ltfat_complex*)c_combined,(const ltfat_complex*)gf_combined,L,W,a,M,
	    f_combined);
   
   mxFree(c_combined);
   mxFree(gf_combined);
   
   plhs[0] = mxCreateDoubleMatrix(L, W, mxCOMPLEX);
   
   combined2split(L*W, (const ltfat_complex*)f_combined, plhs[0]);
   
   mxFree(f_combined);

   return;
   
}


