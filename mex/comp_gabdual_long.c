#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_gabdual_long(g,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   // UGLY, will go away as soon as ltfat backend is unified to a naming convenion
   // Other option is to use forwarder functions
   #undef LTFAT_NAME
   #ifdef LTFAT_SINGLE
   #  define LTFAT_NAME(name) LTFAT_NAME_SINGLE(name)
   #else
   #  define LTFAT_NAME(name) LTFAT_NAME_DOUBLE(name)
   #endif

   int L, R, a, M;

   // Get matrix dimensions.

   L=(int)mxGetM(prhs[0]);
   R=(int)mxGetN(prhs[0]);
   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);

  plhs[0] = ltfatCreateMatrix(L, R,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);
  LTFAT_TYPE* gd_combined = (LTFAT_TYPE*) mxGetData(plhs[0]);
  const LTFAT_TYPE* g_combined = (const LTFAT_TYPE*) mxGetData(prhs[0]);

  #ifdef LTFAT_COMPLEXTYPE
  LTFAT_NAME(gabdual_long)((const LTFAT_REAL (*)[2])g_combined, L, R, a, M, (LTFAT_REAL (*)[2])gd_combined);
  #else
  LTFAT_NAME(gabdualreal_long)(g_combined, L, R, a, M, gd_combined);
  #endif


  return;
}
#endif
/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

// Calling convention:
//  comp_gabdual_long(g,a,M);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, R, a, M;
   ltfat_complex *g_combined, *gd_combined;
   double *gd_r, *gd_i;

   int ii;

   // Get matrix dimensions.

   L=(int)mxGetM(prhs[0]);
   R=(int)mxGetN(prhs[0]);
   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);

   if (mxIsComplex(prhs[0]))
   {
      g_combined  = mxCalloc(L*R,2*sizeof(double));
      gd_combined = mxCalloc(L*R,2*sizeof(double));

      split2combined(L*R, prhs[0], g_combined);

      gabdual_long_d((const ltfat_complex*)g_combined,L,R,a,M,gd_combined);

      plhs[0] = mxCreateDoubleMatrix(L, R, mxCOMPLEX);

      combined2split(L*R, (const ltfat_complex*)gd_combined, plhs[0]);

      mxFree(gd_combined);
      mxFree(g_combined);

   }
   else
   {

      plhs[0] = mxCreateDoubleMatrix(L, R, mxREAL);

      gabdualreal_long_d((const double*)mxGetPr(prhs[0]),
		       L,R,a,M,
		       mxGetPr(plhs[0]));

   }

   return;

}
*/

