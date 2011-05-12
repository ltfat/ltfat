#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_dgtreal_ola(f,g,a,M,bl);
 */


void dothedouble(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[],
		 int L, int W, int a, int M, int N)
{   


}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, gl, W, a, M, N, bl, M2;
   mwSize ndim;
   mwSize dims[3];

   double *f, *g;
   ltfat_complex *out_combined;
   
   /* Get matrix dimensions.*/
   L = mxGetM(prhs[0]); 
   W = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]); 

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   bl=(int)mxGetScalar(prhs[4]);
   
   N=L/a;
   M2=M/2+1;

   dims[0]=M2;
   dims[1]=N;

   if (W==1)
   {
      ndim=2;
   }
   else
   {
      ndim=3;
      dims[2]=W;
   }

   f=mxGetPr(prhs[0]);
   g=mxGetPr(prhs[1]);

   out_combined=(ltfat_complex*)mxMalloc(M*N*W*sizeof(ltfat_complex));

   dgtreal_ola(f,g,L,gl,W,a,M,bl,out_combined);
      
   plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
      
   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);
   return;
   
}


