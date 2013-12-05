#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT

#endif /* _LTFAT_MEX_FILE */ 

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_dwiltiii(f,g,M);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int M, N, L, gl, W;

   // Get matrix dimensions.
   M=(int)mxGetScalar(prhs[2]);
   L=(int)mxGetM(prhs[0]);
   gl=(int) mxGetM(prhs[1]);
   W = mxGetN(prhs[0]);

   N=L/M;

   mwSize dims[]={M, N, W};
   mwSize ndim = W>1?3:2;
   plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);

   const LTFAT_TYPE* f = (const LTFAT_TYPE*) mxGetData(prhs[0]);
   const LTFAT_TYPE* g = (const LTFAT_TYPE*) mxGetData(prhs[1]);
   LTFAT_TYPE* cout = (LTFAT_TYPE*) mxGetData(plhs[0]);

   if(gl<L)
   {
      LTFAT_NAME(dwiltiii_fb)(f,g,L,gl, W, M, cout);
   }
   else
   {
      LTFAT_NAME(dwiltiii_long)(f,g,L,W,M,cout);
   }
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */
