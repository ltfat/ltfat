#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 10
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define REALARGS

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE comp_multidgtrealmp.c
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat/types.h"

// Calling convention:  0 1 2 3     4       5     6     7     8                 9  
//  comp_multidgtrealmp(f,g,a,M,ptype,kernthr,errdb,maxit,maxat,do_pedanticsearch);
//
//  
//

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],
                              int UNUSED(nrhs), const mxArray *prhs[] )
{
    LTFAT_NAME(dgtrealmp_parbuf)* LTFAT_NAME(pbuf) = NULL;
    LTFAT_NAME(dgtrealmp_state)*  LTFAT_NAME(plan) = NULL;
    size_t atoms = 0;
    size_t iters = 0;

    mwSize L  = mxGetNumberOfElements(prhs[0]);
    mwSize dictno = mxGetNumberOfElements(prhs[1]);
    double* aDouble = mxGetData(prhs[2]);
    double* MDouble = mxGetData(prhs[3]);
    int ptype = (int)mxGetScalar(prhs[4]) == 1 ? LTFAT_TIMEINV: LTFAT_FREQINV;;
    double kernthr = mxGetScalar(prhs[5]);
    double errdb = mxGetScalar(prhs[6]);
    size_t maxit = (size_t)mxGetScalar(prhs[7]);
    size_t maxat = (size_t)mxGetScalar(prhs[8]);
    int do_pedanticsearch = (int)mxGetScalar(prhs[9]);

    plhs[0] = mxCreateCellMatrix(dictno, 1);
    LTFAT_COMPLEX** cPtrs = mxMalloc(dictno*sizeof*cPtrs);

    for(int dIdx=0;dIdx<dictno;dIdx++)
    {
        mxSetCell(plhs[0], dIdx,
                  ltfatCreateMatrix(
                      ((mwSize) MDouble[dIdx])/2 + 1, (mwSize)(L/aDouble[dIdx]),
                      LTFAT_MX_CLASSID,mxCOMPLEX));
        cPtrs[dIdx] = mxGetData(mxGetCell(plhs[0],dIdx));
    }

    CHSTAT(LTFAT_NAME(dgtrealmp_parbuf_init)(&LTFAT_NAME(pbuf)));

    for(int dIdx=0;dIdx<dictno;dIdx++)
    {
        CHSTAT(LTFAT_NAME(dgtrealmp_parbuf_add_genwin)(LTFAT_NAME(pbuf),
                mxGetData(mxGetCell(prhs[1],dIdx)),
                mxGetNumberOfElements(mxGetCell(prhs[1],dIdx)),
                (ltfat_int)aDouble[dIdx], (ltfat_int)MDouble[dIdx]));
    }

    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_phaseconv)(LTFAT_NAME(pbuf), ptype));
    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_pedanticsearch)(LTFAT_NAME(pbuf), do_pedanticsearch));
    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_snrdb)(LTFAT_NAME(pbuf), -errdb));
    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_kernrelthr)(LTFAT_NAME(pbuf), kernthr));
    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_maxatoms)(LTFAT_NAME(pbuf), maxat));
    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_maxit)(LTFAT_NAME(pbuf), maxit));
    CHSTAT(LTFAT_NAME(dgtrealmp_setparbuf_iterstep)(LTFAT_NAME(pbuf), L));

    CHSTAT(LTFAT_NAME(dgtrealmp_init)( LTFAT_NAME(pbuf), L, &LTFAT_NAME(plan)));
    CHSTAT(LTFAT_NAME(dgtrealmp_execute_decompose)(LTFAT_NAME(plan), mxGetData(prhs[0]), cPtrs));

    CHSTAT(LTFAT_NAME(dgtrealmp_get_numatoms)(LTFAT_NAME(plan), &atoms));
    CHSTAT(LTFAT_NAME(dgtrealmp_get_numiters)(LTFAT_NAME(plan), &iters));

error:
    if(nlhs>1) plhs[1] = mxCreateDoubleScalar((double)atoms);
    if(nlhs>2) plhs[2] = mxCreateDoubleScalar((double)iters);

    if(LTFAT_NAME(pbuf)) LTFAT_NAME(dgtrealmp_parbuf_done)(&LTFAT_NAME(pbuf));
    if(LTFAT_NAME(plan)) LTFAT_NAME(dgtrealmp_done)(&LTFAT_NAME(plan));
}
#endif /* LTFAT_SINGLE or LTFAT_DOUBLE */

