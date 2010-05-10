#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_dwilt_long(f,g,M,L);
 */


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int M, N, L, W, ii;

   ltfat_complex *f_combined, *g_combined, *cout_combined;
   mwSize ndim;
   mwSize dims[3];

   /* Get matrix dimensions.*/
   M=(int)mxGetScalar(prhs[2]);
   L=(int)mxGetScalar(prhs[3]);
   W = mxGetN(prhs[0]);
    
   N=L/M;

   dims[0]=2*M;
   dims[1]=N/2;

   if (W==1)
   {
      ndim=2;
   }
   else
   {
      ndim=3;
      dims[2]=W;
   }


   if (mxIsComplex(prhs[0]))
   {
      
      f_combined = (ltfat_complex*)mxMalloc(L*W*sizeof(ltfat_complex));
      g_combined = (ltfat_complex*)mxMalloc(L*sizeof(ltfat_complex));
      
      split2combined(L*W, prhs[0], f_combined);

      split2combined(L, prhs[1], g_combined);      
      
      cout_combined = (ltfat_complex*)mxMalloc(L*W*sizeof(ltfat_complex));
      
      dwilt_long((const ltfat_complex*)f_combined,(const ltfat_complex*)g_combined,
		L, W, M, cout_combined);
      mxFree(f_combined);
      mxFree(g_combined);      
      
      plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);

      combined2split(L*W, (const ltfat_complex*)cout_combined,
		     plhs[0]);
      
      mxFree(cout_combined);      
      
   }
   else
   {
      plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
      
      dwiltreal_long(mxGetPr(prhs[0]),mxGetPr(prhs[1]),
		    L, W, M, mxGetPr(plhs[0]));
            
   }

   return;
   
}
