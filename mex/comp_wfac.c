#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_wfac(g,a,M);

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                              int UNUSED(nrhs), const mxArray *prhs[] )
{
   int L, R, N, c, d, p, q;
   ltfatInt a,M,h_a,h_m;

   // Get matrix dimensions.
   L = mxGetM(prhs[0]);
   R = mxGetN(prhs[0]);

   a=(ltfatInt)mxGetScalar(prhs[1]);
   M=(ltfatInt)mxGetScalar(prhs[2]);

   N=L/a;

   c=gcd(a, M, &h_a, &h_m);
   p=a/c;
   q=M/c;
   d=N/q;

   plhs[0] = ltfatCreateMatrix(p*q*R, c*d,LTFAT_MX_CLASSID,mxCOMPLEX);
   LTFAT_COMPLEX* gf_combined = mxGetData(plhs[0]);
   const LTFAT_TYPE* g_combined = mxGetData(prhs[0]);
   LTFAT_NAME(wfac)(g_combined, L, R, a, M, gf_combined);
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */
