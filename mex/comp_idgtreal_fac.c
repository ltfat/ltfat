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
// comp_idgtreal_fac(coef,gf,L,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, W, a, M, N;

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetN(prhs[0])/N;

   plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID,mxREAL);
   const LTFAT_REAL _Complex* c_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   const LTFAT_REAL _Complex* gf_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[1]);
   LTFAT_REAL* f_r = (LTFAT_REAL*) mxGetData(plhs[0]);

   LTFAT_NAME(idgtreal_fac)((const LTFAT_REAL (*)[2])c_combined,
                            (const LTFAT_REAL (*)[2])gf_combined,
                            L,W,a,M,f_r);


   return;
}
#endif

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"


// Calling convention:
// comp_idgtreal_fac(coef,gf,L,a,M);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, W, a, M, N, M2;
   ltfat_complex *gf_combined, *c_combined;
   double *f_r,*gf_r,*gf_i,*c_r,*c_i;

   int ii;

   // This is a floor operation.

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetN(prhs[0])/N;

   M2= M/2+1;

   // Create temporary matrices to convert to correct complex layout.
   c_combined=mxMalloc(M2*N*W*sizeof(ltfat_complex));
   gf_combined=mxMalloc(L*sizeof(ltfat_complex));

   split2combined(M2*N*W, prhs[0], c_combined);
   split2combined(L, prhs[1], gf_combined);

   plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);

   f_r=mxGetPr(plhs[0]);

   LTFAT_NAME_DOUBLE(idgtreal_fac)((const ltfat_complex*)c_combined,
		(const ltfat_complex*)gf_combined,L,W,a,M,
		(double*)f_r);

   mxFree(c_combined);
   mxFree(gf_combined);

   return;

}
*/

