#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 4
#define TYPEDEPARGS 0, 1, 2
#define SINGLEARGS
#define REALARGS

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  cout=comp_gabreassign(s,itime,ifreq,a);

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                              int UNUSED(nrhs), const mxArray *prhs[] )
{
   mwSignedIndex a, M, N, L;
   const LTFAT_REAL *s,*tgrad, *fgrad;
   LTFAT_REAL *sr;

   // Get matrix dimensions.
   M = mxGetM(prhs[0]);
   N = mxGetN(prhs[0]);
   a = (mwSignedIndex)mxGetScalar(prhs[3]);
   L = N*a;

   s     =  mxGetData(prhs[0]);
   tgrad =  mxGetData(prhs[1]);
   fgrad =  mxGetData(prhs[2]);

   plhs[0] = ltfatCreateMatrix(M,N, LTFAT_MX_CLASSID, mxREAL);
   sr      = mxGetData(plhs[0]);

   LTFAT_NAME(gabreassign)(s,tgrad,fgrad,L,1,a,M,sr);
}
#endif
