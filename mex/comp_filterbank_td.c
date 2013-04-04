#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE
/** ISNARGINEQ, ISNARGINLE, ISNARGINGE */
#define ISNARGINEQ 5
/** By defaut, only double version is created. Specifying SINGLEARGS creates second version of the funcntion. */
#define SINGLEARGS 0, 1
/** Specify whether to change the complex number storage format from split planes (Matlab) to interleaved (fftw, complex.h) */
#define TOCOMPLEXFMT 0, 1

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
  #define MEX_FILE __BASE_FILE__
#endif


#include "ltfat_mex_template_helper.h"
#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "mex.h"
#include "math.h"
#include "ltfat.h"
#include "ltfat_types.h"

/*
%COMP_FILTERBANK_TD   Non-uniform filterbank by conv2
%   Usage:  c=comp_filterbank_td(f,g,a,skip,ext);
%
%   Input parameters:
%         f   : Input data - L*W array.
%         g   : Filterbank filters - length M cell-array of vectors of lengths filtLen(m).
%         a   : Subsampling factors - array of length M.
%         skip: Delay of the filters - scalar or array of length M.
%         ext : Border exension technique.
%
%   Output parameters:
%         c  : Cell array of length M. Each element is N(m)*W array.
*/

void TEMPLATE(MEX_FUNC,LTFAT_REAL)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  const mxArray* mxf = prhs[0];
  const mxArray* mxg = prhs[1];
  double* a = (double*) mxGetData(prhs[2]);
  double* skip = (double*) mxGetData(prhs[3]);
  char* ext = mxArrayToString(prhs[4]);

  // input data length
  mwSize L = mxGetM(mxf);
  // number of channels
  mwSize W = mxGetN(mxf);
  // filter number
  mwSize M = mxGetNumberOfElements(mxg);
  // filter lengths
  mwSize filtLen[M];
  //mwSize* filtLen = mxMalloc(M*sizeof(mwSize));
  for(unsigned int m=0;m<M;m++)
  {
     filtLen[m] = (mwSize) mxGetNumberOfElements(mxGetCell(mxg,m));
  }

  // output lengths
  mwSize outLen[M];
  //mwSize* outLen = mxMalloc(M*sizeof(mwSize));
  if(!strcmp(ext,"per"))
  {
     for(unsigned int m = 0; m < M; m++)
     {
        outLen[m] = (mwSize) ceil( L/a[m] );
     }
  }
  else
  {
     for(unsigned int m = 0; m < M; m++)
     {
        outLen[m] = (mwSize) ceil( (L + filtLen[m] - 1 - skip[m] )/a[m] );
     }
  }

  if(mxIsComplex(mxf) || mxIsComplex(mxGetCell(mxg,0)))
  {
     mexErrMsgTxt("Complex inputs not supported yet in a MEX function.");
   //  mxComplexity outComplFlag = mxCOMPLEX;
  }
  else
  {

     mxComplexity outComplFlag = mxREAL;
     //mxClassID outClassID = mxDOUBLE_CLASS;

     // POINTER TO THE INPUT
     LTFAT_REAL* fPtr = (LTFAT_REAL*) mxGetData(prhs[0]);

     // POINTER TO THE FILTERS
     LTFAT_REAL* gPtrs[M];
     // LTFAT_REAL** gPtrs = (LTFAT_REAL**) mxMalloc(M*sizeof(LTFAT_REAL*));
     for(mwIndex m=0;m<M;m++)
     {
        gPtrs[m] = (LTFAT_REAL*) mxGetPr(mxGetCell(mxg, m));
     }

     // POINTER TO OUTPUTS
     LTFAT_REAL* cPtrs[M]; // C99 feature
     //LTFAT_REAL** cPtrs = (LTFAT_REAL**) mxMalloc(M*sizeof(LTFAT_REAL*));
     plhs[0] = mxCreateCellMatrix(M, 1);
     for(mwIndex m=0;m<M;++m)
     {
        mwSize ndim = 2;
        mwSize dims[] = { outLen[m], W};
        /*mwSize dims[2];
        dims[0] =  outLen[m];
        dims[1] =  W;*/
        mxSetCell(plhs[0], m, mxCreateNumericArray(ndim,dims,LTFAT_MX_CLASSID,outComplFlag));
        cPtrs[m] = (LTFAT_REAL*) mxGetData(mxGetCell(plhs[0],m));
        memset(cPtrs[m],0,outLen[m]*W*sizeof(LTFAT_REAL));
     }

     // over all channels
   //  #pragma omp parallel for private(m)

        for(mwIndex m =0; m<M; m++)
        {
          for(mwIndex w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_REAL *fPtrCol = fPtr + w*L;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_REAL *cPtrCol = cPtrs[m] + w*outLen[m];
           //conv_td_sub(fPtrCol,L,&cPtrCol,outLen[m],(const double**)&gPtrs[m],filtLen[m],1,a[m],skip[m],ltfatExtStringToEnum(ext),0);
           LTFAT_NAME(convsub_td)(fPtrCol,L,cPtrCol,outLen[m],gPtrs[m],filtLen[m],a[m],skip[m],ltfatExtStringToEnum(ext));
          }
        }
  }
}
#endif

