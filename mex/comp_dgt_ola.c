#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 6
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

void
LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                         int UNUSED(nrhs), const mxArray *prhs[] )
{
   mwSignedIndex L, gl,W, a, M, N, bl;
   mwSize ndim;
   mwSize dims[3];

   // Get matrix dimensions.
   L = mxGetM(prhs[0]);
   W = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]);

   a=(mwSignedIndex) mxGetScalar(prhs[2]);
   M=(mwSignedIndex) mxGetScalar(prhs[3]);
   bl=(mwSignedIndex) mxGetScalar(prhs[4]);
   mwSignedIndex ptype=(mwSignedIndex) mxGetScalar(prhs[5]);

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
   const LTFAT_COMPLEX* f_combined = mxGetData(prhs[0]);
   const LTFAT_COMPLEX* g_combined = mxGetData(prhs[1]);
   LTFAT_COMPLEX* out_combined = mxGetData(plhs[0]);

   LTFAT_NAME(dgt_ola)(f_combined,g_combined,L,gl,W,a,M,bl,ptype,out_combined);

   return;
}
#endif
