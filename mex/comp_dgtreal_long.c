#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_dgt_long(f,g,a,M);
 */


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, W, a, M, N, ii, M2;
   
   double *f, *g;
   ltfat_complex *out_combined;
   
   /* Get matrix dimensions.*/
   L = mxGetM(prhs[0]); 
   W = mxGetN(prhs[0]);
 
   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   
   N=L/a;
   M2 = (M/2)+1;

   f=mxGetPr(prhs[0]);
   g=mxGetPr(prhs[1]);

   out_combined=(ltfat_complex*)mxMalloc(M2*N*W*sizeof(ltfat_complex));   

   dgtreal_long((const double*)f,(const double*)g,
		L, W, a, M, out_combined);
   
   plhs[0] = mxCreateDoubleMatrix(M2, N*W, mxCOMPLEX);
      
   combined2split(M2*N*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);
   return;
   
}


