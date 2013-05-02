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
//  comp_idgt_fb(coef,g,L,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, W, a, M, N, gl;
   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;

   gl = mxGetM(prhs[1]);
   W  = mxGetM(prhs[0])*mxGetN(prhs[0])/(M*N);

   plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID, mxCOMPLEX);
   const LTFAT_REAL _Complex* c_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   const LTFAT_REAL _Complex* g_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[1]);
   LTFAT_REAL _Complex* f_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

   //#ifdef LTFAT_COMPLEXTYPE
   LTFAT_NAME(idgt_fb)((const LTFAT_REAL (*)[2])c_combined,
                       (const LTFAT_REAL (*)[2])g_combined,
                       L,gl,W,a,M, (LTFAT_REAL (*)[2]) f_combined);

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

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"


// Calling convention:
//  comp_idgt_fb(coef,g,L,a,M);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, W, a, M, N, gl;
   ltfat_complex *f_combined, *g_combined, *c_combined;

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;

   gl = mxGetM(prhs[1]);
   W  = mxGetM(prhs[0])*mxGetN(prhs[0])/(M*N);

   // Create temporary matrices to convert to correct complex layout.

   c_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));
   f_combined=mxMalloc(L*W*sizeof(ltfat_complex));

   split2combined(M*N*W, prhs[0], c_combined);

   if (mxIsComplex(prhs[1]))
   {

      g_combined=mxMalloc(gl*sizeof(ltfat_complex));
      split2combined(   gl, prhs[1], g_combined);

      idgt_fb((const ltfat_complex*)c_combined,(const ltfat_complex*)g_combined,L,gl,W,a,M,
	      f_combined);

      mxFree(g_combined);
   }
   else
   {
      idgt_fb_r((const ltfat_complex*)c_combined,
		(const double*)mxGetPr(prhs[1]),L,gl,W,a,M,
		f_combined);
   }
   mxFree(c_combined);


   plhs[0] = mxCreateDoubleMatrix(L, W, mxCOMPLEX);

   combined2split(L*W, (const ltfat_complex*)f_combined, plhs[0]);

   mxFree(f_combined);

   return;

}
*/

