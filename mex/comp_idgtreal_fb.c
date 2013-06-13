#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5

#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXARGS

#define MATCHEDARGS 1

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_idgtreal_fb(coef,g,L,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, W, a, M, N, gl;

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetN(prhs[0])/N;
   gl = mxGetM(prhs[1]);


   plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID,mxREAL);
   const LTFAT_REAL _Complex* c_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   const LTFAT_REAL * gf = (const LTFAT_REAL *) mxGetData(prhs[1]);
   LTFAT_REAL* f_r = (LTFAT_REAL*) mxGetData(plhs[0]);

   LTFAT_NAME(idgtreal_fb)((const LTFAT_REAL (*)[2])c_combined,
                            gf,L,gl,W,a,M,f_r);


   return;
}
#endif

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"


// Calling convention:
//  comp_idgtreal_fb(coef,g,L,a,M);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, W, a, M, N, M2, gl;
   ltfat_complex *c_combined;

   int ii;

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetN(prhs[0])/N;
   gl = mxGetM(prhs[1]);


   M2= M/2+1;

   // Create temporary matrices to convert to correct complex layout.
   c_combined=mxMalloc(M2*N*W*sizeof(ltfat_complex));

   split2combined(M2*N*W, prhs[0], c_combined);

   plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);

   LTFAT_NAME_DOUBLE(idgtreal_fb)((const ltfat_complex*)c_combined,
		(const double*)mxGetPr(prhs[1]),L,gl,W,a,M,
		(double*)mxGetPr(plhs[0]));

   mxFree(c_combined);

   return;

}
*/

