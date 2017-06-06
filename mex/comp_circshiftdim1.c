#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINGE 1
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define NOCOMPLEXFMTCHANGE

#endif /* _LTFAT_MEX_FILE */

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
#define MEX_FILE __BASE_FILE__
#endif
#include "ltfat_mex_template_helper.h"


#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

/*
COMP_CELLCOEF2TF Cell to a tf-layout
   Usage: coef = comp_cellcoef2tf(coef,maxLen)
*/

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray* plhs[],
                              int nrhs, const mxArray* prhs[] )
{
    const mxArray* mxin = prhs[0];
    LTFAT_REAL* inPtr = mxGetData(mxin);
    mwSignedIndex L = mxGetM(mxin);
    mwSize W = mxGetN(mxin);
    mwSignedIndex shift = (mwSignedIndex) mxGetScalar(prhs[1]);

    if (mxIsComplex(mxin))
        plhs[0] = ltfatCreateNdimArray(mxGetNumberOfDimensions(mxin),
                                       mxGetDimensions(mxin),
                                       LTFAT_MX_CLASSID, mxCOMPLEX);
    else
        plhs[0] = ltfatCreateNdimArray(mxGetNumberOfDimensions(mxin),
                                       mxGetDimensions(mxin),
                                       LTFAT_MX_CLASSID, mxREAL);


    LTFAT_REAL* outPtr = mxGetData(plhs[0]);

    mwSignedIndex p = (L - shift) % L;

    if (p < 0) p += L;

    for (mwSize w = 0; w < W; w++)
    {
        const LTFAT_REAL* inPtrCol = inPtr + w * L;
        LTFAT_REAL* outPtrCol = outPtr + w * L;

        memcpy(outPtrCol, inPtrCol + p, (L - p)*sizeof * outPtrCol);
        memcpy(outPtrCol + L - p, inPtrCol, p * sizeof * outPtrCol);
    }

    if (mxIsComplex(mxin))
    {
        LTFAT_REAL* inImagPtr = mxGetImagData(mxin);
        LTFAT_REAL* outImagPtr = mxGetImagData(plhs[0]);
        for (mwSize w = 0; w < W; w++)
        {
            const LTFAT_REAL* inPtrCol = inImagPtr + w * L;
            LTFAT_REAL* outPtrCol = outImagPtr + w * L;

            memcpy(outPtrCol, inPtrCol + p, (L - p)*sizeof * outPtrCol);
            memcpy(outPtrCol + L - p, inPtrCol, p * sizeof * outPtrCol);
        }

    }


}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */
