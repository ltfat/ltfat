#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_dgtreal_fb(f,g,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, gl,W, a, M, N, M2;
   ltfat_complex *out_combined;
   double *f, *g;
   
   int ii;
   
   /* Get matrix dimensions.*/
   L  = mxGetM(prhs[0]); 
   W  = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]); 

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   M2=M/2+1;

   N=L/a;
   
   /* Create temporary matrices to convert to correct complex layout. */
   
   out_combined=mxMalloc(M2*N*W*sizeof(ltfat_complex));
   
   f=mxGetPr(prhs[0]);
   g=mxGetPr(prhs[1]);
   
   dgtreal_fb((const double*)f,(const double*)g,L,gl,W,a,M,
              out_combined);
      
   plhs[0] = mxCreateDoubleMatrix(M2, N*W, mxCOMPLEX);
   
   combined2split(M2*N*W, (const ltfat_complex*)out_combined, plhs[0]);

      
   mxFree(out_combined);
   return;
   
}


