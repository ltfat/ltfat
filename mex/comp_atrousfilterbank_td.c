#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define NARGINEQ 4
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
#include "ltfat_types.h"

/*
%COMP_ATROUSFILTERBANK_TD   Uniform filterbank by conv2
%   Usage:  c=comp_atrousfilterbank_fft(f,g,a,offset);
%
%   Input parameters:
%         f   : Input data - L*W array.
%         g   : Filterbank filters - filtLen*M array.
%         a   : Filter upsampling factor - scalar.
%         offset: Delay of the filters - scalar or array of length M.
%
%   Output parameters:
%         c  : L*M*W array of coefficients
%
*/
void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    const mxArray* mxf = prhs[0];
    const mxArray* mxg = prhs[1];
    double* a = mxGetPr(prhs[2]);
    double* offset = mxGetPr(prhs[3]);

    // input data length
    unsigned int L = mxGetM(mxf);
    // number of channels
    unsigned int W = mxGetN(mxf);
    // filter number
    unsigned int M = mxGetN(mxg);
    // filter length
    unsigned int filtLen = mxGetM(mxg);


        // POINTER TO THE INPUT
        LTFAT_TYPE* fPtr = (LTFAT_TYPE*) mxGetPr(prhs[0]);

        // POINTER TO THE FILTERS
        LTFAT_TYPE** gPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
        for(unsigned int m=0; m<M; m++)
        {
            gPtrs[m] = ((LTFAT_TYPE*) mxGetData(mxg)) + m*filtLen;
        }

        mwSize ndim = 3;
        mwSize dims[3];
        dims[0] =  L;
        dims[1] =  M;
        dims[2] =  W;
        plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);

        // over all filters
        for(unsigned int m =0; m<M; m++)
        {
            // over all channels
            for(unsigned int w =0; w<W; w++)
            {
                // Obtain pointer to w-th column in input
                LTFAT_TYPE *fPtrCol = fPtr + w*L;
                LTFAT_TYPE *cPtrPlane = ((LTFAT_TYPE*) mxGetData(plhs[0])) + w*L*M;
                // Obtaing pointer to w-th column in m-th element of output cell-array
                LTFAT_TYPE *cPtrCol = cPtrPlane + m*L;
                LTFAT_NAME(atrousconvsub_td)(fPtrCol, L, cPtrCol, L, gPtrs[m], filtLen, (int)*a, -offset[m], ltfatExtStringToEnum("per"));
            }
        }
}
#endif


