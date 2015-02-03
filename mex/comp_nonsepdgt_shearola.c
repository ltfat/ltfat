#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 8
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_dgt_ola(f,g,a,M,s0,s1,br,bl);

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                              int UNUSED(nrhs), const mxArray *prhs[] )
{

   mwSize ndim;
   mwSize dims[3];

   // Get matrix dimensions.
   const int L = mxGetM(prhs[0]);
   const int W = mxGetN(prhs[0]);
   const int gl = mxGetM(prhs[1]);

   const int a=(int)mxGetScalar(prhs[2]);
   const int M=(int)mxGetScalar(prhs[3]);

   const int s0=(int)mxGetScalar(prhs[4]);
   const int s1=(int)mxGetScalar(prhs[5]);
   const int br=(int)mxGetScalar(prhs[6]);

   const int bl=(int)mxGetScalar(prhs[7]);

   const int N=L/a;

   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,mxCOMPLEX);
   const LTFAT_COMPLEX* f_combined = (LTFAT_COMPLEX*) mxGetData(prhs[0]);
   const LTFAT_COMPLEX* g_combined = (LTFAT_COMPLEX*) mxGetData(prhs[1]);
   LTFAT_COMPLEX* out_combined = (LTFAT_COMPLEX*) mxGetData(plhs[0]);

   LTFAT_NAME(dgt_shearola)(f_combined,g_combined,L,gl,W,a,M,s0,s1,br,bl,out_combined);

   return;
}
#endif

