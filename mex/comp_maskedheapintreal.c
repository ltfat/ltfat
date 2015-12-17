#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 9
#define TYPEDEPARGS 0
#define MATCHEDARGS 1, 2, 8
#define SINGLEARGS
#define COMPLEXARGS

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
// phase=comp_maskedheapintreal(s,itime,ifreq,mask,a,M,tol,do_timeinv,usephase);

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
        int UNUSED(nrhs), const mxArray *prhs[] )
{
    int a, M, N, L, W, M2, phasetype, useoutphase;
    double tol;

    const LTFAT_COMPLEX* c;
    const LTFAT_REAL *tgrad, *fgrad;
    const double* maskDouble;
    LTFAT_REAL *phase;
    int* mask;

    // Get inputs
    c     = mxGetData(prhs[0]);
    tgrad = mxGetData(prhs[1]);
    fgrad = mxGetData(prhs[2]);
    maskDouble  = mxGetData(prhs[3]);
    a     = (int)mxGetScalar(prhs[4]);
    tol   = mxGetScalar(prhs[6]);
    phasetype = (int)mxGetScalar(prhs[7]);
    mxArray* knownphase = prhs[8];

    useoutphase = mxGetNumberOfElements(knownphase) > 0;

    // Get matrix dimensions.
    M = (int) mxGetScalar(prhs[5]);
    M2 = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    L = N * a;
    W = 1;

    mask = ltfat_malloc(M2 * N * W * sizeof * mask);

    for (ltfatInt w = 0; w < M2 * N * W; w++ )
        mask[w] = (int) maskDouble[w];

    // Create output matrix and zero it.
    plhs[0] = ltfatCreateMatrix(M2, N, LTFAT_MX_CLASSID, mxREAL);

    // Get pointer to output
    phase = mxGetData(plhs[0]);

    if(useoutphase)
        memcpy(phase,mxGetData(knownphase),M2 * N * W * sizeof * phase);

    LTFAT_NAME(maskedheapintreal)(c, tgrad, fgrad, mask, a, M, L, W, tol,
                                  phasetype, useoutphase, phase );

    ltfat_free(mask);
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */

