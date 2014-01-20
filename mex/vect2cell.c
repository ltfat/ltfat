#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define NOCOMPLEXFMTCHANGE

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Assuming __BASE_FILE__ is known by the compiler.
   Otherwise specify this filename
   e.g. #define MEX_FILE "comp_col2diag.c"  */
#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// LTFAT_COMPLEXTYPE, LTFAT_SINGLE, LTFAT_DOUBLE

/* Calling convention:
 *  c=vect2cell(x,idx);
 */
void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    mwSize L = mxGetM(prhs[0]);
    mwSize W = mxGetN(prhs[0]);
    mwSize M = mxGetNumberOfElements(prhs[1]);
    double* Lc = mxGetData(prhs[1]);

    plhs[0] = mxCreateCellMatrix(M,1);
    LTFAT_REAL* cPr[M];
    LTFAT_REAL* cPi[M];

    for(mwIndex ii=0;ii<M;ii++)
    {
        mxArray* tmpA = ltfatCreateMatrix((mwSize)Lc[ii], W,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY);
        mxSetCell(plhs[0],ii,tmpA);
        cPr[ii] = (LTFAT_REAL*) mxGetPr(tmpA);
        cPi[ii] = (LTFAT_REAL*) mxGetPi(tmpA);
    }

    LTFAT_REAL* xPr = (LTFAT_REAL*) mxGetPr(prhs[0]);
    LTFAT_REAL* xPi = (LTFAT_REAL*) mxGetPi(prhs[0]);


    for(mwIndex w=0;w<W;w++)
    {
       LTFAT_REAL* xTmp = xPr + w*L;
       for(mwIndex ii=0;ii<M;ii++)
       {
           mwSize LcTmp = (mwSize)Lc[ii];
           LTFAT_REAL* cTmp = cPr[ii] + w*LcTmp;
           memcpy(cTmp,xTmp,LcTmp*sizeof(LTFAT_REAL));
           xTmp+=LcTmp;
       }
    }

      #ifdef LTFAT_COMPLEXTYPE
    for(mwIndex w=0;w<W;w++)
    {
       LTFAT_REAL* xTmp = xPi + w*L;
       for(mwIndex ii=0;ii<M;ii++)
       {
           mwSize LcTmp = (mwSize)Lc[ii];
           LTFAT_REAL* cTmp = cPi[ii] + w*LcTmp;
           memcpy(cTmp,xTmp,LcTmp*sizeof(LTFAT_REAL));
           xTmp+=LcTmp;
       }
    }
      #endif

}
#endif

