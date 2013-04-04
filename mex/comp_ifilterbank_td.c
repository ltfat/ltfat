#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

/*  Define which arguments are to be checked and cast to single if either of them is single. */
#define NARGINEQ 6
#define DATATYPECHECK 0, 1

/* Specify whether to change the complex number storage format from split planes (Matlab) to interleaved (fftw, complex.h) */
//#define CHCOMPLEXFORMAT 1

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
  #define MEX_FILE __BASE_FILE__
//#else
//#define MEX_FILE "comp_ifilterbank_td.c"
#endif

/* The following header includes this file twice setting either LTFAT_SINGLE or LTFAT_DOUBLE.
    At the end of the header, LTFAT_SINGLE or LTFAT_DOUBLE is unset. */
#include "ltfat_mex_template_helper.h"
/* Do not allow processing this file further unless LTFAT_SINGLE or LTFAT_DOUBLE is specified. */
#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)


/** From now on, it is like ordinary mexFunction but params prhs[PRHSTOCHECK[i]], i=0:length(PRHSTOCHECK) are now of type LTFAT_REAL.
    Complex array still has to be processed separatelly.
    Enclose calls to the ltfat backend in LTFAT_NAME() macro.
    Enclose calls to the fftw in LTFAT_FFTW() macro.
    Avoid using mx functions working with concrete data type.
    e.g. use mxGetData intead of mxGetPr (or recast to LTFAT_REAL*)
         mxCreateNumericArray with macro LTFAT_MX_CLASSID instead of createDoubleMatrix
 */
#include "mex.h"
#include "math.h"
#include "ltfat.h"
/**  The following defines single and double versions for the types and macros:
  LTFAT_COMPLEX - fftw_complex or fftwf_complex
  LTFAT_REAL - double or float
  LTFAT_NAME(name) - for LTFAT_SINGLE add "s" to the beginning of the function name
  LTFAT_FFTW(name) - adds "fftw_" or "fftwf_" to the beginning of the function name
  LTFAT_MX_CLASSID - mxDOUBLE_CLASS or mxSINGLE_CLASS
**/
#include "ltfat_types.h"

/*
%COMP_IFILTERBANK_TD   Synthesis filterbank by conv2
%   Usage:  f=comp_ifilterbank_fft(c,g,a,Ls,skip,ext);
%
%   Input parameters:
%         c    : Cell array of length M, each element is N(m)*W matrix.
%         g    : Filterbank filters - length M cell-array, each element is vector of length filtLen(m)
%         a    : Upsampling factors - array of length M.
%         skip : Delay of the filters - scalar or array of length M.
%         Ls   : Output length.
%         ext  : Border exension technique.
%
%   Output parameters:
%         f  : Output Ls*W array.
%
*/
void TEMPLATE(MEX_FUNC,LTFAT_REAL)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  // printf("Filename: %s, Function name %s, %d \n.",__FILE__,__func__,mxIsDouble(prhs[0]));
  const mxArray* mxc = prhs[0];
  const mxArray* mxg = prhs[1];
  double* a = mxGetPr(prhs[2]);
  double* Lsdouble = mxGetPr(prhs[3]);
  unsigned int Ls = (unsigned int) *Lsdouble;
  double* skip = mxGetPr(prhs[4]);
  char* ext = mxArrayToString(prhs[5]);

  // number of channels
  unsigned int W = mxGetN(mxGetCell(mxc,0));

  // filter number
  unsigned int M = mxGetNumberOfElements(mxg);

  // input data length
  unsigned int* Lc = mxMalloc(M*sizeof(unsigned int));
  for(unsigned int m=0;m<M;m++)
  {
     Lc[m] = (unsigned int) mxGetM(mxGetCell(mxc,m));
  }

  // filter lengths
  unsigned int* filtLen = mxMalloc(M*sizeof(unsigned int));
  for(unsigned int m=0;m<M;m++)
  {
     filtLen[m] = (unsigned int) mxGetNumberOfElements(mxGetCell(mxg,m));
  }

  if(mxIsComplex(mxc) || mxIsComplex(mxGetCell(mxg,0)))
  {
     mexErrMsgTxt("Complex inputs not supported yet in a MEX function.");
   //  mxComplexity outComplFlag = mxCOMPLEX;
  }
  else
  {
     mxComplexity outComplFlag = mxREAL;
     //mxClassID outClassID = mxDOUBLE_CLASS;

     // POINTER TO THE INPUT
     LTFAT_REAL** cPtrs = (LTFAT_REAL**) mxMalloc(M*sizeof(LTFAT_REAL*));
     for(unsigned int m=0;m<M;++m)
     {
        cPtrs[m] = (LTFAT_REAL*) mxGetData(mxGetCell(mxc,m));
     }

     // allocate output
     mwSize ndim = 2;
     mwSize dims[2];
     dims[0] =  Ls;
     dims[1] =  W;
     plhs[0] = mxCreateNumericArray(ndim,dims,LTFAT_MX_CLASSID,outComplFlag);


      // POINTER TO OUTPUT
     LTFAT_REAL* fPtr = (LTFAT_REAL*) mxGetData(plhs[0]);
     // Set to zeros
     memset(fPtr,0,Ls*W*sizeof(LTFAT_REAL));

     // POINTER TO THE FILTERS
     LTFAT_REAL** gPtrs = (LTFAT_REAL**) mxMalloc(M*sizeof(LTFAT_REAL*));
     for(unsigned int m=0;m<M;m++)
     {
        gPtrs[m] = (LTFAT_REAL*) mxGetData(mxGetCell(mxg, m));
     }

     // over all channels
   //  #pragma omp parallel for private(m)

        for(unsigned int m =0; m<M; m++)
        {
          for(unsigned int w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_REAL *fPtrCol = fPtr + w*Ls;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_REAL *cPtrCol = cPtrs[m] + w*Lc[m];
           //(upconv_td)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts, int fLen, int up, int skip, enum ltfatWavExtType ext)
           LTFAT_NAME(upconv_td)(cPtrCol,Lc[m],fPtrCol,Ls,gPtrs[m],filtLen[m],a[m],skip[m],ltfatExtStringToEnum(ext));
          }
       }
  }
}
#endif
