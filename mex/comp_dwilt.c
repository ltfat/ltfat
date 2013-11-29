#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_dwilt(f,g,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int M, N, L, gl, W;

   mwSize ndim;
   mwSize dims[3];

   // Get matrix dimensions.
   M=(int)mxGetScalar(prhs[2]);
   L=(int)mxGetM(prhs[0]);
   gl=(int) mxGetM(prhs[1]);
   W = mxGetN(prhs[0]);

   N=L/M;

   dims[0]=2*M;
   dims[1]=N/2;

   if (W==1)
   {
      ndim=2;
   }
   else
   {
      ndim=3;
      dims[2]=W;
   }

   plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);

   const LTFAT_TYPE* f = (const LTFAT_TYPE*) mxGetData(prhs[0]);
   const LTFAT_TYPE* g = (const LTFAT_TYPE*) mxGetData(prhs[1]);
   LTFAT_TYPE* cout = (LTFAT_TYPE*) mxGetData(plhs[0]);

   if(gl<L)
   {
      LTFAT_NAME(dwilt_fb)(f,g,L,gl, W, M, cout);
   }
   else
   {
      LTFAT_NAME(dwilt_long)(f,g,L,W,M,cout);
   }
}
#endif
