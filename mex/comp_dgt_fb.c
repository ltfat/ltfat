#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 4
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_dgt_fb(f,g,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L  = mxGetM(prhs[0]);
   int W  = mxGetN(prhs[0]);
   int gl = mxGetM(prhs[1]);

   int a=(int)mxGetScalar(prhs[2]);
   int M=(int)mxGetScalar(prhs[3]);

   int N=L/a;

   mwSize ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   mwSize dims[3] = {M, N, W};
   plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,mxCOMPLEX);
   const LTFAT_TYPE* f_combined = (const LTFAT_TYPE*) mxGetData(prhs[0]);
   const LTFAT_TYPE* g_combined = (const LTFAT_TYPE*) mxGetData(prhs[1]);
   LTFAT_COMPLEX* out_combined = (LTFAT_COMPLEX*) mxGetData(plhs[0]);

   LTFAT_NAME(dgt_fb)(f_combined, g_combined,L,gl,W,a,M,out_combined);
}
#endif


/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

// Calling convention:
// comp_dgt_fb(f,g,a,M);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, gl,W, a, M, N;
   mwSize ndim;
   mwSize dims[3];

   ltfat_complex *f_combined, *g_combined;
   ltfat_complex *out_combined;

   // Get matrix dimensions.
   L  = mxGetM(prhs[0]);
   W  = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]);

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);

   N=L/a;

   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   out_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));

   // Copy the data.

   if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
   {
      f_combined=mxMalloc(L*W*sizeof(ltfat_complex));
      split2combined(L*W, prhs[0], f_combined);

      g_combined=mxMalloc(L*2*sizeof(ltfat_complex));
      split2combined(gl,  prhs[1], g_combined);

      dgt_fb((const ltfat_complex*)f_combined,(const ltfat_complex*)g_combined,
	     L,gl,W,a,M,
	     (ltfat_complex*)out_combined);

      mxFree(f_combined);
      mxFree(g_combined);

   }
   else
   {
      dgt_fb_r((const double*)mxGetPr(prhs[0]),
	       (const double*)mxGetPr(prhs[1]),
	       L,gl,W,a,M,
	       (ltfat_complex*)out_combined);
   }

   plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);

   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);

   return;

}

*/
