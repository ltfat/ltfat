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
   int L, W, gl, a, M, N;
   mwSize ndim;
   mwSize dims[3];

   ltfat_complex *f_combined, *g_combined, *out_combined;
   
   /* Get matrix dimensions.*/
   L = mxGetM(prhs[0]); 
   W = mxGetN(prhs[0]);

   gl = mxGetM(prhs[1]); 
   M  = mxGetN(prhs[1]);
 
   a=(int)mxGetScalar(prhs[2]);
   
   N=L/a;

   dims[0]=N;
   dims[1]=M;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   f_combined=(ltfat_complex*)mxMalloc(L*W*sizeof(ltfat_complex));
   g_combined=mxMalloc(gl*M*sizeof(ltfat_complex));
   out_combined=(ltfat_complex*)mxMalloc(N*M*W*sizeof(ltfat_complex));
   
   split2combined(L*W,  prhs[0], f_combined);
   split2combined(gl*M, prhs[1], g_combined);

   ufilterbank_fft((const ltfat_complex*)f_combined,(const ltfat_complex*)g_combined,
		   L,gl,W,a,M,out_combined);
   
   mxFree(f_combined);   
   mxFree(g_combined);
   
   plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
      
   combined2split(N*M*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);
   return;
   
}


