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
//  comp_wfac(g,a,M);

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
   int L, R, a, M, N, c, d, p, q,h_a,h_m;

   // Get matrix dimensions.
   L = mxGetM(prhs[0]);
   R = mxGetN(prhs[0]);

   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);

   N=L/a;

   c=gcd(a, M,&h_a, &h_m);
   p=a/c;
   q=M/c;
   d=N/q;

  plhs[0] = ltfatCreateMatrix(p*q*R, c*d,LTFAT_MX_CLASSID,mxCOMPLEX);
  LTFAT_REAL _Complex* gf_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

  #ifdef LTFAT_COMPLEXTYPE
  const LTFAT_REAL _Complex* g_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
  LTFAT_NAME(wfac)((const LTFAT_REAL (*)[2])g_combined, L, R, a, M, (LTFAT_REAL (*)[2])gf_combined);
  #else
  const LTFAT_REAL* g_combined = (const LTFAT_REAL*) mxGetData(prhs[0]);
  LTFAT_NAME(wfac_r)(g_combined, L, R, a, M, (LTFAT_REAL (*)[2])gf_combined);
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
//  comp_wfac(g,a,M);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, R, a, M, N, c, d, p, q,h_a,h_m;
   ltfat_complex *g_combined, *gf_combined;

   // Get matrix dimensions.
   L = mxGetM(prhs[0]);
   R = mxGetN(prhs[0]);

   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);

   N=L/a;

   c=gcd(a, M,&h_a, &h_m);
   p=a/c;
   q=M/c;
   d=N/q;

   // Create temporary matrices to convert to correct complex layout.
   gf_combined=mxMalloc(L*R*sizeof(ltfat_complex));

   if (mxIsComplex(prhs[0]))
   {
      // g is complex
      g_combined=mxMalloc(L*R*sizeof(ltfat_complex));

      split2combined(L*R, prhs[0], g_combined);

      wfac((const ltfat_complex*)g_combined,
	   L, R, a, M, (ltfat_complex*)gf_combined);

      mxFree(g_combined);
   }
   else
   {
      // g is real
     wfac_r(mxGetPr(prhs[0]), L, R, a, M,
	    gf_combined);
   }

   plhs[0] = mxCreateDoubleMatrix(p*q*R, c*d, mxCOMPLEX);

   combined2split(L*R, (const ltfat_complex*)gf_combined, plhs[0]);

   mxFree(gf_combined);

   return;

}
*/

