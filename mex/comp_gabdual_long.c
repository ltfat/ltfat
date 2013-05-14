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

   int L, R, a, M;

   // Get matrix dimensions.

   L=(int)mxGetM(prhs[0]);
   R=(int)mxGetN(prhs[0]);
   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);

  plhs[0] = ltfatCreateMatrix(L, R,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);
  LTFAT_TYPE* gd_combined = (LTFAT_TYPE*) mxGetData(plhs[0]);
  const LTFAT_TYPE* g_combined = (const LTFAT_TYPE*) mxGetData(prhs[0]);

  LTFAT_NAME(gabdual_long)(g_combined, L, R, a, M, gd_combined);

  return;
}
#endif
