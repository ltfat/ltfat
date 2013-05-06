#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 4
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define REALARGS

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
// comp_dgtreal_fb(f,g,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, gl,W, a, M, N, M2;

   // Get matrix dimensions.
   L  = mxGetM(prhs[0]);
   W  = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]);

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   M2=M/2+1;

   N=L/a;

   plhs[0] = ltfatCreateMatrix(M2, N*W,LTFAT_MX_CLASSID,mxCOMPLEX);
   const LTFAT_REAL * f = (const LTFAT_REAL *) mxGetData(prhs[0]);
   const LTFAT_REAL * g = (const LTFAT_REAL *) mxGetData(prhs[1]);
   LTFAT_REAL _Complex* out_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

   LTFAT_NAME(dgtreal_fb)(f,g,L,gl,W,a,M,(LTFAT_REAL (*)[2])out_combined);

   return;
}
#endif

/*

#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

// Calling convention:
// comp_dgtreal_fb(f,g,a,M);


 void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, gl,W, a, M, N, M2;
   ltfat_complex *out_combined;
   double *f, *g;

   // Get matrix dimensions.
   L  = mxGetM(prhs[0]);
   W  = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]);

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   M2=M/2+1;

   N=L/a;

   // Create temporary matrices to convert to correct complex layout.

   out_combined=mxMalloc(M2*N*W*sizeof(ltfat_complex));

   f=mxGetPr(prhs[0]);
   g=mxGetPr(prhs[1]);

   dgtreal_fb((const double*)f,(const double*)g,L,gl,W,a,M,
              out_combined);

   plhs[0] = mxCreateDoubleMatrix(M2, N*W, mxCOMPLEX);

   combined2split(M2*N*W, (const ltfat_complex*)out_combined, plhs[0]);


   mxFree(out_combined);
   return;

}

*/
