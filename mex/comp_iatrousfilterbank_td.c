#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

/*  Define which arguments are to be checked and cast to single if either of them is single. */
static const int PRHSTOCHECK[] = { 0, 1};

/* Specify whether to change the complex number storage format from split planes (Matlab) to interleaved (fftw, complex.h) */
//#define CHCOMPLEXFORMAT 1

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
#define MEX_FILE __BASE_FILE__
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
%COMP_IATROUSFILTERBANK_TD   Synthesis Uniform filterbank by conv2
%   Usage:  f=comp_iatrousfilterbank_fft(c,g,a,skip);
%
%   Input parameters:
%         c    : L*M*W array of coefficients.
%         g    : Filterbank filters - filtLen*M array.
%         a    : Filters upsampling factor - scalar.
%         skip : Delay of the filters - scalar or array of length M.
%
%   Output parameters:
%         f  : Output L*W array.
%
*/
void TEMPLATE(MEX_FUNC, LTFAT_REAL)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    // printf("Filename: %s, Function name %s, %d \n.",__FILE__,__func__,mxIsDouble(prhs[0]));
    const mxArray* mxc = prhs[0];
    const mxArray* mxg = prhs[1];
    double* a = mxGetPr(prhs[2]);
    double* skip = mxGetPr(prhs[3]);

    // number of channels
    const mwSize *dims = mxGetDimensions(mxc);
    unsigned int L = dims[0];
    unsigned int M = dims[1];
    unsigned int W = 1;
    if(mxGetNumberOfDimensions(mxc)>2)
    {
        W = dims[2];
    }

    // filter number
    //unsigned int M = mxGetNumberOfElements(mxg);

    // filter length
    unsigned int filtLen = mxGetM(mxg);

    if(mxIsComplex(mxc) || mxIsComplex(mxg))
    {
        mexErrMsgTxt("Complex inputs not supported yet in a MEX function.");
        //  mxComplexity outComplFlag = mxCOMPLEX;
    }
    else
    {
        mxComplexity outComplFlag = mxREAL;
        //mxClassID outClassID = mxDOUBLE_CLASS;

        // allocate output
        mwSize ndim = 2;
        mwSize dims[2];
        dims[0] =  L;
        dims[1] =  W;
        plhs[0] = mxCreateNumericArray(ndim,dims,LTFAT_MX_CLASSID,outComplFlag);

        // POINTER TO OUTPUT
        LTFAT_REAL* fPtr = (LTFAT_REAL*) mxGetData(plhs[0]);
        // Set to zeros
        memset(fPtr,0,L*W*sizeof(LTFAT_REAL));

        // POINTER TO THE FILTERS
        LTFAT_REAL** gPtrs = (LTFAT_REAL**) mxMalloc(M*sizeof(LTFAT_REAL*));
        for(unsigned int m=0; m<M; m++)
        {
            gPtrs[m] = ((LTFAT_REAL*) mxGetData(mxg)) + m*filtLen;
        }

        // over all channels
        //  #pragma omp parallel for private(m)
        for(unsigned int m =0; m<M; m++)
        {
            for(unsigned int w =0; w<W; w++)
            {
                // Obtain pointer to w-th column in input
                LTFAT_REAL *fPtrCol = fPtr + w*L;
                LTFAT_REAL *cPtrPlane = ((LTFAT_REAL*) mxGetData(mxc)) + w*L*M;

                // Obtaing pointer to w-th column in m-th element of output cell-array
                LTFAT_REAL *cPtrCol = cPtrPlane + m*L;
                //(upconv_td)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts, int fLen, int up, int skip, enum ltfatWavExtType ext)
                LTFAT_NAME(atrousupconv_td)(cPtrCol,L,fPtrCol,L,gPtrs[m],filtLen,(int)*a,skip[m],ltfatExtStringToEnum("per"));
            }
        }
    }
}
#endif
