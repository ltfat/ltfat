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
//  comp_dgt_ola(f,g,a,M,bl);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int L, gl,W, a, M, N, bl;
   mwSize ndim;
   mwSize dims[3];

   // Get matrix dimensions.
   L = mxGetM(prhs[0]);
   W = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]);

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   bl=(int)mxGetScalar(prhs[4]);

   N=L/a;

   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,mxCOMPLEX);
   const LTFAT_REAL _Complex* f_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   const LTFAT_REAL _Complex* g_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[1]);
   LTFAT_REAL _Complex* out_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

   LTFAT_NAME(dgt_ola)((const LTFAT_REAL (*)[2])f_combined,
                       (const LTFAT_REAL (*)[2])g_combined,
                        L,gl,W,a,M,bl,(LTFAT_REAL (*)[2])out_combined);

   return;
}
#endif
