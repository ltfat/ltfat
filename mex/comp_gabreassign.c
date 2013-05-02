#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 4
#define TYPEDEPARGS 0, 1, 2
#define SINGLEARGS
#define REALARGS

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  cout=comp_gabreassign(s,itime,ifreq,a);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int a, M, N, L;
   const LTFAT_REAL *s,*tgrad, *fgrad;
   LTFAT_REAL *sr;

   // Get matrix dimensions.
   M = mxGetM(prhs[0]);
   N = mxGetN(prhs[0]);
   a = (int)mxGetScalar(prhs[3]);
   L = N*a;

   s     = (const LTFAT_REAL*) mxGetPr(prhs[0]);
   tgrad = (const LTFAT_REAL*) mxGetPr(prhs[1]);
   fgrad = (const LTFAT_REAL*) mxGetPr(prhs[2]);

   plhs[0] = ltfatCreateMatrix(M,N, LTFAT_MX_CLASSID, mxREAL);
   sr      = (LTFAT_REAL*) mxGetPr(plhs[0]);

   LTFAT_NAME(gabreassign)(s,tgrad,fgrad,L,1,a,M,sr);

   return;
}
#endif

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"

// Calling convention:
//  cout=comp_gabreassign(s,itime,ifreq,a);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int a, M, N, L;
   double *s,*tgrad, *fgrad,*sr;

   // Get matrix dimensions.
   M = mxGetM(prhs[0]);
   N = mxGetN(prhs[0]);
   a = (int)mxGetScalar(prhs[3]);
   L = N*a;

   s     = mxGetPr(prhs[0]);
   tgrad = mxGetPr(prhs[1]);
   fgrad = mxGetPr(prhs[2]);

   plhs[0] = mxCreateDoubleMatrix(M,N, mxREAL);
   sr      = mxGetPr(plhs[0]);

   gabreassign((double*)s,(double*)tgrad,(double*)fgrad,
               L,1,a,M,(double*)sr);

   return;

}
*/

