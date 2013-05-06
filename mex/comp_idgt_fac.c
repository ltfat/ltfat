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
//  comp_idgt_fac(coef,gf,L,a,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, W, a, M, N;

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetM(prhs[0])*mxGetN(prhs[0])/(M*N);


   plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID, mxCOMPLEX);
   const LTFAT_REAL _Complex* c_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   const LTFAT_REAL _Complex* gf_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[1]);
   LTFAT_REAL _Complex* f_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

   LTFAT_NAME(idgt_fac)((const LTFAT_REAL (*)[2])c_combined,
                        (const LTFAT_REAL (*)[2])gf_combined,
                        L,W,a,M,(LTFAT_REAL (*)[2])f_combined);



   return;
}
#endif

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

// Calling convention:
//  comp_idgt_fac(coef,gf,L,a,M);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int L, W, a, M, N;
   ltfat_complex *f_combined, *gf_combined, *c_combined;

   // Get matrix dimensions.
   L=(int)mxGetScalar(prhs[2]);
   a=(int)mxGetScalar(prhs[3]);
   M=(int)mxGetScalar(prhs[4]);
   N=L/a;
   W = mxGetM(prhs[0])*mxGetN(prhs[0])/(M*N);

   // Create temporary matrices to convert to correct complex layout.

   c_combined=mxMalloc(M*N*W*sizeof(ltfat_complex));
   gf_combined=mxMalloc(L*sizeof(ltfat_complex));
   f_combined=mxMalloc(L*W*sizeof(ltfat_complex));

   split2combined(M*N*W, prhs[0], c_combined);
   split2combined(L, prhs[1], gf_combined);

   idgt_fac((const ltfat_complex*)c_combined,(const ltfat_complex*)gf_combined,L,W,a,M,
	    f_combined);

   mxFree(c_combined);
   mxFree(gf_combined);

   plhs[0] = mxCreateDoubleMatrix(L, W, mxCOMPLEX);

   combined2split(L*W, (const ltfat_complex*)f_combined, plhs[0]);

   mxFree(f_combined);

   return;

}
*/

