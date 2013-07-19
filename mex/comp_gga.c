#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
//#define NOCOMPLEXFMTCHANGE

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)



// Calling convention:
//  c = comp_gga(f,indvec)

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

   mwSize L  = mxGetM(prhs[0]);
   mwSize W  = mxGetN(prhs[0]);
   mwSize M = mxGetNumberOfElements(prhs[1]);

   const LTFAT_TYPE* fPtr = (const LTFAT_TYPE*) mxGetPr(prhs[0]);
   const double* indVecPtr = (const double*) mxGetPr(prhs[1]);

   plhs[0] = ltfatCreateMatrix(M,W,LTFAT_MX_CLASSID,mxCOMPLEX);
   LTFAT_REAL _Complex* cPtr = (LTFAT_REAL _Complex*) mxGetPr(plhs[0]);

   LTFAT_NAME(gga)(fPtr,indVecPtr,L,W,M,cPtr);


   return;
}
#endif