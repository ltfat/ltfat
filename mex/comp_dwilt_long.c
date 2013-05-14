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
//  comp_dwilt_long(f,g,M,L);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   // UGLY, but there are no cd_dwilt_long cs_dwilt_long functions
   // Other option is to use forwarder functions
   #undef LTFAT_NAME
   #ifdef LTFAT_SINGLE
   #  define LTFAT_NAME(name) LTFAT_NAME_SINGLE(name)
   #else
   #  define LTFAT_NAME(name) LTFAT_NAME_DOUBLE(name)
   #endif

   int M, N, L, W;

   mwSize ndim;
   mwSize dims[3];

   // Get matrix dimensions.
   M=(int)mxGetScalar(prhs[2]);
   L=(int)mxGetScalar(prhs[3]);
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


   #ifdef LTFAT_COMPLEXTYPE
      LTFAT_NAME(dwilt_long)((const LTFAT_REAL (*)[2])f,
                             (const LTFAT_REAL (*)[2])g,
                              L, W, M,
                             (LTFAT_REAL (*)[2]) cout);
   #else
      LTFAT_NAME(dwiltreal_long)(f,g,L, W, M, cout);
   #endif


   return;
}
#endif
