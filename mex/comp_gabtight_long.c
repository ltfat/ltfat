#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_gabtight_long(g,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int L, a, M;
   ltfat_complex *g_combined, *gd_combined;
   double *gd_r, *gd_i;
   
   int ii;
   
   /* Get matrix dimensions.*/
   
   L=(int)mxGetM(prhs[0]);
   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);

   if (mxIsComplex(prhs[0]))
   {
      g_combined  = mxMalloc(L*sizeof(ltfat_complex));
      gd_combined = mxMalloc(L*sizeof(ltfat_complex));

      split2combined(L, prhs[0], g_combined);

      gabtight_long((const ltfat_complex*)g_combined,L,a,M,gd_combined);

      plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);
   
      combined2split(L, (const ltfat_complex*)gd_combined, plhs[0]);

      mxFree(gd_combined);
      mxFree(g_combined);

   }
   else      
   {

      plhs[0] = mxCreateDoubleMatrix(L, 1, mxREAL);

      gabtightreal_long((const double*)mxGetPr(prhs[0]),
		       L,a,M,
		       mxGetPr(plhs[0]));
   
   }

   return;
  
}


