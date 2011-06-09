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
   mwSize ndim;
   mwSize dims[3];

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

   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }
   
   out_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));
   
   /* Copy the data. */
         
   if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
   {
      f_combined=mxMalloc(L*W*sizeof(ltfat_complex));
      split2combined(L*W, prhs[0], f_combined);

      g_combined=mxMalloc(L*2*sizeof(ltfat_complex));
      split2combined(gl,  prhs[1], g_combined);

      dgt_fb((const ltfat_complex*)f_combined,(const ltfat_complex*)g_combined,
	     L,gl,W,a,M,
	     (ltfat_complex*)out_combined);

      mxFree(f_combined);
      mxFree(g_combined);
      
   }
   else
   {
      dgt_fb_r((const double*)mxGetPr(prhs[0]),
	       (const double*)mxGetPr(prhs[1]),
	       L,gl,W,a,M,
	       (ltfat_complex*)out_combined);
   }

   plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);   

   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);
   
   mxFree(out_combined);
   
   return;
   
}


