#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 8
#define TYPEDEPARGS 0, 1, 2, 3
#define SINGLEARGS
#define REALARGS

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//                            0 1     2     3     4 5 6   7
// phase=comp_heapintreal_ufb(s,tgrad,fgrad,cfreq,a,M,tol,phasetype);
// phasetype defines how to adjust tgrad and fgrad such that
// phase corresponds to:
// phasetype  0:  freqinv
// phasetype  1:  timeinv
// phasetype  2:  do not adjust the gradient, it is already correct

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray* plhs[],
                              int UNUSED(nrhs), const mxArray* prhs[] )
{
    // Get inputs
    const mxArray* mxs  = prhs[0];
    const LTFAT_REAL* s = mxGetData(mxs);
    const LTFAT_REAL* tgrad = mxGetData(prhs[1]);
    const LTFAT_REAL* fgrad = mxGetData(prhs[2]);
    const LTFAT_REAL* cfreq = mxGetData(prhs[3]);
    mwSize a  = (mwSize)mxGetScalar(prhs[4]);
    mwSize M = (mwSize)mxGetScalar(prhs[5]);
    LTFAT_REAL tol = (LTFAT_REAL) mxGetScalar(prhs[6]);
    int phasetype  = (int)mxGetScalar(prhs[7]);

    // Get matrix dimensions.
    mwSize N = ltfatGetN(prhs[0]);
    mwSize L = N * a;
    mwSize W = 1;
    phasetype--;

    if (mxGetNumberOfDimensions(mxs) > 2)
        W = mxGetDimensions(mxs)[2];

    // Create output matrix and zero it.
    plhs[0] = ltfatCreateNdimArray(mxGetNumberOfDimensions(mxs),
                                   mxGetDimensions(mxs),
                                   LTFAT_MX_CLASSID, mxREAL);

    // Get pointer to output
    LTFAT_REAL* phase = mxGetData(plhs[0]);

    if (phasetype == 1)
        LTFAT_NAME(heapintreal_ufb)(s, tgrad, fgrad, cfreq, a, M, L, W, tol, phase);
    else
        LTFAT_NAME(heapintreal_relgrad_ufb)(s, tgrad, fgrad, cfreq, a, M, L, W, tol,
                                        phase);
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */
