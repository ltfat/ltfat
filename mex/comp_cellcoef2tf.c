#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
  #define MEX_FILE __BASE_FILE__
#endif
#include "ltfat_mex_template_helper.h"


#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
/** USER DEFINED HEADERS **/
#include "math.h"

/*
COMP_CELLCOEF2TF Cell to a tf-layout
   Usage: coef = comp_cellcoef2tf(coef,dim)
*/

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

  const mxArray* mxCoef = prhs[0];
  //const mxArray* mxDim = prhs[1];


  mwSize M = mxGetNumberOfElements(mxCoef);
  
  mwSize maxCoefElLen = 0;
  mwSize coefElLen[M];
  const LTFAT_TYPE* coefElPtr[M];
  LTFAT_TYPE* coefOutPtr;
  
  for(mwIndex ii=0;ii<M;ii++)
  {
     mxArray* mxCoefEl = mxGetCell(mxCoef,ii); 
	 coefElPtr[ii] = (const LTFAT_TYPE*) mxGetPr(mxCoefEl);
	 coefElLen[ii] = mxGetM(mxCoefEl);
	 if(maxCoefElLen<coefElLen[ii])
	 {
	    maxCoefElLen = coefElLen[ii];
	 }
  }

  plhs[0] = ltfatCreateMatrix(M, maxCoefElLen ,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);
  coefOutPtr = (LTFAT_TYPE*) mxGetPr(plhs[0]);
  
  for(mwIndex m=0;m<M;m++)
  {
     const LTFAT_TYPE* coefElPtrTmp = coefElPtr[m];
	 LTFAT_TYPE* coefOutPtrTmp = coefOutPtr+m;
	 double lenRatio = ((double)coefElLen[m]-1)/((double)maxCoefElLen-1);
	 for(mwIndex ii=0;ii<maxCoefElLen;ii++)
	 {
	    *coefOutPtrTmp = coefElPtrTmp[(mwIndex)((ii*lenRatio)+0.5)];
		coefOutPtrTmp += M;
	 }
  }
}
#endif

