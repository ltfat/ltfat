#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 6
#define TYPEDEPARGS 0, 1, 2
#define SINGLEARGS
#define REALARGS

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
// phase=comp_heapint(s,itime,ifreq,a,tol);

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                              int UNUSED(nrhs), const mxArray *prhs[] )
{
    int a, M, N, L, W, M2;
    double tol;

    const LTFAT_REAL *s, *tgrad, *fgrad;
    LTFAT_REAL *phase;

    // Get inputs
    s     = mxGetData(prhs[0]);
    tgrad = mxGetData(prhs[1]);
    fgrad = mxGetData(prhs[2]);
    a     = (int)mxGetScalar(prhs[3]);
    tol   = mxGetScalar(prhs[5]);

    // Get matrix dimensions.
    M = (int)mxGetScalar(prhs[4]);
    M2 = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    L = N * a;
    W = 1;

    // Create output matrix and zero it.
    plhs[0] = ltfatCreateMatrix(M2, N, LTFAT_MX_CLASSID, mxREAL);

    // Get pointer to output
    phase = (LTFAT_REAL*) mxGetPr(plhs[0]);

    LTFAT_NAME(heapintreal)(s, tgrad, fgrad, a, M, L, W, tol, phase ); 
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */

