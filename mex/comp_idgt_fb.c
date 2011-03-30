#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"


/* Calling convention:
 *  comp_idgt_fb(coef,g,L,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, W, a, M, N, gl;
   ltfat_complex *f_combined, *g_combined, *c_combined;   
   
   /* Get matrix dimensions.*/
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;

   gl = mxGetM(prhs[1]); 
   W  = mxGetM(prhs[0])*mxGetN(prhs[0])/(M*N); 

   /* Create temporary matrices to convert to correct complex layout. */

   c_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));   
   f_combined=mxMalloc(L*W*sizeof(ltfat_complex));
   
   split2combined(M*N*W, prhs[0], c_combined);

   if (mxIsComplex(prhs[1]))
   {

      g_combined=mxMalloc(gl*sizeof(ltfat_complex));
      split2combined(   gl, prhs[1], g_combined);
      
      idgt_fb((const ltfat_complex*)c_combined,(const ltfat_complex*)g_combined,L,gl,W,a,M,
	      f_combined);
      
      mxFree(g_combined);     
   }
   else
   {
      idgt_fb_r((const ltfat_complex*)c_combined,
		(const double*)mxGetPr(prhs[1]),L,gl,W,a,M,
		f_combined);            
   }
   mxFree(c_combined);

   
   plhs[0] = mxCreateDoubleMatrix(L, W, mxCOMPLEX);
   
   combined2split(L*W, (const ltfat_complex*)f_combined, plhs[0]);
   
   mxFree(f_combined);

   return;
   
}


