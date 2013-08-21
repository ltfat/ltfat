#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_filterbank_fft

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "math.h"

// Calling convention:
// c = comp_filterbank_fft(F,G,a)

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  const mxArray* mxF = prhs[0];
  const mxArray* mxG = prhs[1];
  double* a = (double*) mxGetData(prhs[2]);

  // input data length
  mwSize L = mxGetM(mxF);
  // number of channels
  mwSize W = mxGetN(mxF);
  // filter number
  mwSize M = mxGetNumberOfElements(mxG);

  // output lengths
  mwSize outLen[M];
  //mwSize* outLen = mxMalloc(M*sizeof(mwSize));

  for(unsigned int m = 0; m < M; m++)
  {
     outLen[m] = (mwSize) ceil( L/a[m] );
  }
 

     // POINTER TO THE INPUT
     LTFAT_REAL _Complex* FPtr = (LTFAT_REAL _Complex*) mxGetData(prhs[0]);

     // POINTER TO THE FILTERS
     LTFAT_REAL _Complex* GPtrs[M];
     // LTFAT_TYPE** gPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     for(mwIndex m=0;m<M;m++)
     {
        GPtrs[m] = (LTFAT_REAL _Complex*) mxGetPr(mxGetCell(mxG, m));
     }

     // POINTER TO OUTPUTS
     LTFAT_REAL _Complex* cPtrs[M]; // C99 feature
     //LTFAT_TYPE** cPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     plhs[0] = mxCreateCellMatrix(M, 1);
     for(mwIndex m=0;m<M;++m)
     {
        mxSetCell(plhs[0], m, ltfatCreateMatrix(outLen[m], W,LTFAT_MX_CLASSID,mxCOMPLEX));
        cPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(plhs[0],m));
        memset(cPtrs[m],0,outLen[m]*W*sizeof(LTFAT_REAL _Complex));
     }

     // over all channels
   //  #pragma omp parallel for private(m)

        for(mwIndex m =0; m<M; m++)
        {
          for(mwIndex w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_REAL _Complex *FPtrCol = FPtr + w*L;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_REAL _Complex *cPtrCol = cPtrs[m] + w*outLen[m];
          
           LTFAT_NAME(convsub_fft)(FPtrCol,GPtrs[m],L,(size_t)a[m],cPtrCol);
          }
        }
}
#endif




