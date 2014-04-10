#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5
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

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],
                              int nrhs, const mxArray *prhs[] )
{
    int a, M, N, L, W;
    double tol;

    const LTFAT_REAL *s, *tgrad, *fgrad;
    LTFAT_REAL *phase;

    // Get inputs
    s     = (const LTFAT_REAL*) mxGetPr(prhs[0]);
    tgrad = (const LTFAT_REAL*) mxGetPr(prhs[1]);
    fgrad = (const LTFAT_REAL*) mxGetPr(prhs[2]);
    a     = (int)mxGetScalar(prhs[3]);
    tol   = mxGetScalar(prhs[4]);

    // Get matrix dimensions.
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    L = N * a;
    W = 1;

    // Create output matrix and zero it.
    plhs[0] = ltfatCreateMatrix(M, N, LTFAT_MX_CLASSID, mxREAL);

    // Get pointer to output
    phase = (LTFAT_REAL*) mxGetPr(plhs[0]);

    LTFAT_NAME(heapint)(s, tgrad, fgrad, a, M, L, W, tol, phase);
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */



/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"

// Calling convention:
// phase=comp_heapint(s,itime,ifreq,a,tol);


void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[] )

{
int a, M, N, L, W;
double tol;

double *s, *tgrad, *fgrad, *phase;

// Get inputs
s     = mxGetPr(prhs[0]);
tgrad = mxGetPr(prhs[1]);
fgrad = mxGetPr(prhs[2]);
a     = (int)mxGetScalar(prhs[3]);
tol   = mxGetScalar(prhs[4]);

// Get matrix dimensions.
M = mxGetM(prhs[0]);
N = mxGetN(prhs[0]);
L = N*a;
W = 1;

// Create output matrix and zero it.
plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);

// Get pointer to output
phase=mxGetPr(plhs[0]);

heapint((const double*)s,
(const double*)tgrad,
(const double*)fgrad,
a, M, L, W,
tol, phase);

return;

}
*/

