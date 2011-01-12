#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_dgt_fb(f,g,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, gl,W, a, M, N;
   ltfat_complex *f_combined, *g_combined;
   ltfat_complex *out_combined;
   
   int ii;
   
   /* Get matrix dimensions.*/
   L  = mxGetM(prhs[0]); 
   W  = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]); 

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   
   N=L/a;
   
   /* Create temporary matrices to convert to correct complex layout. */
   
   f_combined=mxMalloc(L*W*sizeof(ltfat_complex));
   g_combined=mxMalloc(L*2*sizeof(ltfat_complex));
   out_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));
   
   /* Copy the data. */
   split2combined(L*W, prhs[0], f_combined);
   split2combined(gl,  prhs[1], g_combined);
         
   dgt_fb((const ltfat_complex*)f_combined,(const ltfat_complex*)g_combined,
	  L,gl,W,a,M,
	  (ltfat_complex*)out_combined);
   
   mxFree(f_combined);
   mxFree(g_combined);
   
   plhs[0] = mxCreateDoubleMatrix(M, N*W, mxCOMPLEX);

   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);
   
   mxFree(out_combined);
   
   return;
   
}


