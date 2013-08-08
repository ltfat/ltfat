#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

/*  Define which arguments are to be checked and cast to single if either of them is single. */
#define NARGINEQ 6
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT

/* Specify whether to change the complex number storage format from split planes (Matlab) to interleaved (fftw, complex.h) */
//#define CHCOMPLEXFORMAT 1

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
  #define MEX_FILE __BASE_FILE__
//#else
//#define MEX_FILE "comp_ifilterbank_td.c"
#endif


#include "ltfat_mex_template_helper.h"
#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)

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
%COMP_IFILTERBANK_TD   Synthesis filterbank
%   Usage:  f=comp_ifilterbank_fft(c,g,a,Ls,offset,ext);
%
%   Input parameters:
%         c    : Cell array of length M, each element is N(m)*W matrix.
%         g    : Filterbank filters - length M cell-array, each element is vector of length filtLen(m)
%         a    : Upsampling factors - array of length M.
%         offset : Delay of the filters - scalar or array of length M.
%         Ls   : Output length.
%         ext  : Border exension technique.
%
%   Output parameters:
%         f  : Output Ls*W array.
%
*/
void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  // printf("Filename: %s, Function name %s, %d \n.",__FILE__,__func__,mxIsDouble(prhs[0]));
  const mxArray* mxc = prhs[0];
  const mxArray* mxg = prhs[1];
  double* a = mxGetPr(prhs[2]);
  double* Lsdouble = mxGetPr(prhs[3]);
  unsigned int Ls = (unsigned int) *Lsdouble;
  double* offset = mxGetPr(prhs[4]);
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

     // POINTER TO THE INPUT
     LTFAT_TYPE** cPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     for(unsigned int m=0;m<M;++m)
     {
        cPtrs[m] = (LTFAT_TYPE*) mxGetData(mxGetCell(mxc,m));
     }

     // allocate output
     plhs[0] = ltfatCreateMatrix(Ls, W,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);


      // POINTER TO OUTPUT
     LTFAT_TYPE* fPtr = (LTFAT_TYPE*) mxGetData(plhs[0]);
     // Set to zeros
     memset(fPtr,0,Ls*W*sizeof(LTFAT_TYPE));

     // POINTER TO THE FILTERS
     LTFAT_TYPE** gPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     for(unsigned int m=0;m<M;m++)
     {
        gPtrs[m] = (LTFAT_TYPE*) mxGetData(mxGetCell(mxg, m));
     }

     // over all channels
   //  #pragma omp parallel for private(m)

        for(unsigned int m =0; m<M; m++)
        {
          for(unsigned int w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_TYPE *fPtrCol = fPtr + w*Ls;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_TYPE *cPtrCol = cPtrs[m] + w*Lc[m];
           //(upconv_td)(const LTFAT_TYPE *in, int inLen, LTFAT_TYPE *out, const int outLen, const LTFAT_TYPE *filts, int fLen, int up, int skip, enum ltfatWavExtType ext)
           LTFAT_NAME(upconv_td)(cPtrCol,Lc[m],fPtrCol,Ls,gPtrs[m],filtLen[m],a[m],-offset[m],ltfatExtStringToEnum(ext));
          }
       }

}
#endif
