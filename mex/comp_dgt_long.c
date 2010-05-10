#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_dgt_long(f,g,a,M);
 */


void dothedouble(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[],
		 int L, int W, int a, int M, int N)
{   


}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
  int L, W, a, M, N, ii;

  double *f_r,*f_i;
   ltfat_complex *gf, *f_combined, *g_combined, *out_combined;
   
   /* Get matrix dimensions.*/
   L = mxGetM(prhs[0]); 
   W = mxGetN(prhs[0]);
 
   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   
   N=L/a;

   f_r=mxGetPr(prhs[0]);
   f_i=mxGetPi(prhs[0]);

   gf = mxMalloc(L*sizeof(ltfat_complex));

   /* Create factorization of window */   
   if (mxIsComplex(prhs[1]))
   {
      /* g is complex */
      g_combined=mxMalloc(L*sizeof(ltfat_complex));
      split2combined(L, prhs[1], g_combined);
      wfac((const ltfat_complex*)g_combined, L, a, M, gf);
      mxFree(g_combined);
   }
   else
   {
      /* g is real */
      wfac_r(mxGetPr(prhs[1]), L, a, M, gf);
   }

   f_combined=(ltfat_complex*)mxMalloc(L*W*sizeof(ltfat_complex));
   out_combined=(ltfat_complex*)mxMalloc(M*N*W*sizeof(ltfat_complex));
   
   split2combined(L*W, prhs[0], f_combined);

   dgt_fac((const ltfat_complex*)f_combined,
	   (const ltfat_complex*)gf,L,W,a,M,out_combined);
   
   mxFree(f_combined);
   mxFree(gf);
   
   plhs[0] = mxCreateDoubleMatrix(M, N*W, mxCOMPLEX);
      
   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);
   return;
   
}


