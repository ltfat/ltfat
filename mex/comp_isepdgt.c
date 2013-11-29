#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_isepdgt(coef,g,L,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, W, a, M, N, gl;
   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;

   gl = mxGetM(prhs[1]);
   W  = mxGetNumberOfElements(prhs[0])/(M*N);

   plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID, mxCOMPLEX);
   const LTFAT_COMPLEX* c_combined = (const LTFAT_COMPLEX*) mxGetData(prhs[0]);
   const LTFAT_COMPLEX* g_combined = (const LTFAT_COMPLEX*) mxGetData(prhs[1]);
   LTFAT_COMPLEX* f_combined = (LTFAT_COMPLEX*) mxGetData(plhs[0]);

   if(gl<L)
   {
      LTFAT_NAME(idgt_fb)(c_combined,g_combined, L,gl,W,a,M, f_combined);
   }
   else
   {
      LTFAT_NAME(idgt_long)(c_combined,g_combined, L,W,a,M, f_combined); 
   }
  /* #else
   NOT CALLING idgt_fb_r:
   TO DO: Do it better.
   LTFAT_NAME(idgt_fb_r)((const LTFAT_REAL (*)[2])c_combined,
                         g_combined,
                         L,gl,W,a,M,(LTFAT_REAL (*)[2]) f_combined);
   #endif
   */
   return;
}
#endif
