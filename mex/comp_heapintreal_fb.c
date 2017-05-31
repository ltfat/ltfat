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
//                           0 1     2     3     4    5     6 7 8 9         10  11 
// phase=comp_heapintreal_fb(s,tgrad,fgrad,neigh,dist,cfreq,a,M,N,chanStart,tol,phasetype);

// PHASETYPE IS NOT USED!
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
    const mxArray* mxneigh  = prhs[3];
    const mxArray* mxdist  = prhs[4];
    
    const LTFAT_REAL* s = mxGetData(mxs);

    const LTFAT_REAL* tgrad = mxGetData(prhs[1]);
    const LTFAT_REAL* fgrad = mxGetData(prhs[2]);  
    const LTFAT_REAL* cfreq = mxGetData(prhs[5]);
    
    LTFAT_REAL* a = mxGetData(prhs[6]);  
    mwSize M = (mwSize)mxGetScalar(prhs[7]);

    const ltfatInt* N  = mxGetData(prhs[8]);
    const ltfatInt* chanStart  = mxGetData(prhs[9]);

    LTFAT_REAL tol = (LTFAT_REAL) mxGetScalar(prhs[10]);
    int phasetype = (int)mxGetScalar(prhs[11]);

    const ltfatInt sLen = chanStart[M+1];
   
    const LTFAT_REAL* neighPtr[W*4*chanStart[sLen]];
    const LTFAT_REAL* distPtr[W*3*chanStart[sLen]];

    // Get matrix dimensions.
    ltfatInt W = 1;

    if (mxGetNumberOfDimensions(mxs) > 2)
        W = mxGetDimensions(mxs)[2];

    for ( int w = 0; w < W; w++ ) 
    {
        for ( int m = 0; m < sLen; m++ ) 
        {
            // Create entries of distPtr
            distPtr[m+(3*w)*sLen] = mxGetData(mxGetCell(mxdist, m+(3*w)*sLen));
            distPtr[m+(3*w+1)*sLen] = mxGetData(mxGetCell(mxdist, m+(3*w+1)*sLen));
            distPtr[m+(3*w+2)*sLen] = mxGetData(mxGetCell(mxdist, m+(3*w+2)*sLen));

            // Create entries of neighPtr
            neighPtr[m+(4*w)*sLen] = mxGetData(mxGetCell(mxdist, m+(4*w)*sLen));
            neighPtr[m+(4*w+1)*sLen] = mxGetData(mxGetCell(mxdist, m+(4*w+1)*sLen));
            neighPtr[m+(4*w+2)*sLen] = mxGetData(mxGetCell(mxdist, m+(4*w+2)*sLen));
            neighPtr[m+(4*w+3)*sLen] = mxGetData(mxGetCell(mxdist, m+(4*w+3)*sLen));
        } 
    }
      
    // Create output matrix and zero it.
    plhs[0] = ltfatCreateNdimArray(mxGetNumberOfDimensions(mxs),
                                   mxGetDimensions(mxs),
                                   LTFAT_MX_CLASSID, mxREAL);

    // Get pointer to output
    LTFAT_REAL* phase = mxGetData(plhs[0]);

    mwSize L = N[0] * a[0];     
  
    if (phasetype == 1)
        LTFAT_NAME(heapint_fb)(s, tgrad, fgrad, neighPtr, distPtr, cfreq, 
			a, M, N, chanStart, L, W, tol, phase);
    else
        LTFAT_NAME(heapint_relgrad_fb)(s, tgrad, fgrad, neighPtr, distPtr,
			cfreq, a, M, N, chanStart, L, W, tol, phase);
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */
