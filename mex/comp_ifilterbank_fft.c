#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_ifilterbank_fft

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "math.h"
#include "config.h"

// Calling convention:
// c = comp_ifilterbank_fft(c,G,a)

static LTFAT_FFTW(plan)* LTFAT_NAME(oldPlans) = 0;
static mwSize* LTFAT_NAME(oldLc) = 0;

static mwSize LTFAT_NAME(oldM) = 0;

void LTFAT_NAME(ifftMexAtExitFnc)()
{
    if(LTFAT_NAME(oldPlans)!=0)
    {
        for(mwIndex m=0; m<LTFAT_NAME(oldM); m++)
        {
            if(LTFAT_NAME(oldPlans)[m]!=0)
            {
                LTFAT_FFTW(destroy_plan)(LTFAT_NAME(oldPlans)[m]);
            }
        }
        ltfat_free(LTFAT_NAME(oldPlans));
        LTFAT_NAME(oldPlans) = 0;
    }

    if(LTFAT_NAME(oldLc)!=0)
    {
        ltfat_free(LTFAT_NAME(oldLc));
        LTFAT_NAME(oldLc) = 0;
    }
    LTFAT_NAME(oldM) = 0;
#ifdef _DEBUG
    mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
#endif
}

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],
                              int nrhs, const mxArray *prhs[] )
{
    static int atExitFncRegistered = 0;
    if(!atExitFncRegistered)
    {
        LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(ifftMexAtExitFnc));
        atExitFncRegistered = 1;
    }

    const mxArray* mxc = prhs[0];
    const mxArray* mxG = prhs[1];
    double* aDouble = (double*) mxGetData(prhs[2]);

    // input data length
    mwSize L = mxGetM(mxGetCell(mxG,0));
    // number of channels
    mwSize W = mxGetN(mxGetCell(mxc,0));
    // filter number
    mwSize M = mxGetNumberOfElements(mxc);

    // input lengths
    mwSize inLen[M];

    // Hop sizes array
    mwSize a[M];

    plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID,mxCOMPLEX);

    // POINTER TO THE OUTPUT
    LTFAT_COMPLEX* FPtr = (LTFAT_COMPLEX*) mxGetData(plhs[0]);


    // POINTERS TO THE FILTERS
    LTFAT_COMPLEX* GPtrs[M];

    // POINTER TO INPUTS
    LTFAT_COMPLEX* cPtrs[M]; // C99 feature
    LTFAT_COMPLEX* cbuf[M]; // C99 feature

    if(M!=LTFAT_NAME(oldM))
    {
        LTFAT_NAME(ifftMexAtExitFnc)();
        LTFAT_NAME(oldM) = M;
        LTFAT_NAME(oldLc) = (mwSize*) ltfat_calloc(M,sizeof(mwSize));
        LTFAT_NAME(oldPlans) = ltfat_calloc(M,sizeof(LTFAT_FFTW(plan)));
    }

    for(mwIndex m =0; m<M; m++)
    {
        a[m] = (mwSize) aDouble[m];
        GPtrs[m] = (LTFAT_COMPLEX*) mxGetData(mxGetCell(mxG, m));
        cPtrs[m] = (LTFAT_COMPLEX*) mxGetData(mxGetCell(mxc,m));
        inLen[m] = mxGetM(mxGetCell(mxc,m));
        cbuf[m] = ltfat_malloc(inLen[m]*sizeof*cbuf[m]);

        if(LTFAT_NAME(oldLc)[m]!=inLen[m])
        {
            LTFAT_NAME(oldLc)[m] = inLen[m];

            LTFAT_FFTW(plan) ptmp = LTFAT_FFTW(plan_dft_1d)(inLen[m], cbuf[m], cbuf[m], FFTW_FORWARD, FFTW_OPTITYPE);

            if(LTFAT_NAME(oldPlans)[m]!=0)
            {
                LTFAT_FFTW(destroy_plan)(LTFAT_NAME(oldPlans)[m]);
            }
            LTFAT_NAME(oldPlans)[m]=ptmp;
        }
    }

    LTFAT_NAME(ifilterbank_fft_plans)(cPtrs,GPtrs, L, W, a, M, FPtr,  LTFAT_NAME(oldPlans), cbuf);

    for(mwIndex m =0; m<M; m++)
    {
        ltfat_free(cbuf[m]);
    }

}
#endif




