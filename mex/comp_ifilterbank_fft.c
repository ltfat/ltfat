#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_ifilterbank_fft

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "math.h"

// Calling convention:
// c = comp_ifilterbank_fft(c,G,a)

void LTFAT_NAME(ifftMexAtExitFnc)()
{
   #ifdef _DEBUG
   mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
   #endif
}

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  #ifdef _DEBUG
  static int atExitFncRegistered = 0;
  if(!atExitFncRegistered)
  {
     LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(ifftMexAtExitFnc));
     atExitFncRegistered = 1;
  }
  #endif

  const mxArray* mxc = prhs[0];
  const mxArray* mxG = prhs[1];
  double* a = (double*) mxGetData(prhs[2]);

  // input data length
  mwSize L = mxGetM(mxGetCell(mxG,0));
  // number of channels
  mwSize W = mxGetN(mxGetCell(mxG,0));
  // filter number
  mwSize M = mxGetNumberOfElements(mxG);

  // input lengths
  mwSize cLen[M];
 
  plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID,mxCOMPLEX);

  // POINTER TO THE OUTPUT
  LTFAT_REAL _Complex* FPtr = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);
  memset(FPtr,0,L*W*sizeof(LTFAT_REAL _Complex));

  // POINTERS TO THE FILTERS
  LTFAT_REAL _Complex* GPtrs[M];
	 
  // POINTER TO INPUTS
  LTFAT_REAL _Complex* cPtrs[M]; // C99 feature

  for(mwIndex m=0;m<M;m++)
  {
     GPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(mxG, m));
     cPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(mxc,m));
	 cLen[m] = mxGetM(mxGetCell(mxc,m));
  }

  // over all channels
  //  #pragma omp parallel for private(m)

  for(mwIndex m =0; m<M; m++)
  {
    for(mwIndex w =0; w<W; w++)
    {
       // Obtain pointer to w-th column in output
       LTFAT_REAL _Complex *FPtrCol = FPtr + w*L;
       // Obtaing pointer to w-th column in m-th element of input cell-array
       LTFAT_REAL _Complex *cPtrCol = cPtrs[m] + w*cLen[m];
          
       LTFAT_NAME(upconv_fft)(cPtrCol,cLen[m],GPtrs[m],(size_t)a[m],FPtrCol);
    }
  }
}
#endif




