#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 6
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define REALARGS

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
//  comp_dgtreal_ola(f,g,a,M,bl);

void
LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                         int UNUSED(nrhs), const mxArray *prhs[] )
{
    mwSignedIndex L, gl, W, a, M, N, bl, M2;
    mwSize ndim;
    mwSize dims[3];

    // Get matrix dimensions.
    L = mxGetM(prhs[0]);
    W = mxGetN(prhs[0]);
    gl = mxGetM(prhs[1]);

    a=(mwSignedIndex)mxGetScalar(prhs[2]);
    M=(mwSignedIndex)mxGetScalar(prhs[3]);
    bl=(mwSignedIndex)mxGetScalar(prhs[4]);
    mwSignedIndex phasetype=(mwSignedIndex)mxGetScalar(prhs[5]);

    N=L/a;
    M2=M/2+1;

    dims[0]=M2;
    dims[1]=N;

    if (W==1)
    {
        ndim=2;
    }
    else
    {
        ndim=3;
        dims[2]=W;
    }

    plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID, mxCOMPLEX);
    const LTFAT_REAL * f = mxGetData(prhs[0]);
    const LTFAT_REAL * g = mxGetData(prhs[1]);
    LTFAT_COMPLEX* out_combined = (LTFAT_COMPLEX*) mxGetData(plhs[0]);

    LTFAT_NAME(dgtreal_ola)(f,g,L,gl,W,a,M,bl,phasetype,out_combined);
    return;
}
#endif
