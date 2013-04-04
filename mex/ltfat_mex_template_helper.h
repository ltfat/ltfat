//#define MEX_FILE
/** Allow including this file only if MEX_FILE is set */
#if defined(MEX_FILE)

/** Allow including this file only once */
#ifndef _LTFAT_MEX_TEMPLATEHELPER_H
#define _LTFAT_MEX_TEMPLATEHELPER_H 1

/** General mex function name. */
#define MEX_FUNC ltfat_mexFunction

/** Template macros */
#define CAT(X,Y) X##_##Y
#define TEMPLATE(X,Y) CAT(X,Y)

/** Undefine LTFAT_DOUBLE and LTFAT_SINGLE if they are set */
#ifdef LTFAT_DOUBLE
#undef LTFAT_DOUBLE
#endif
#ifdef LTFAT_SINGLE
#undef LTFAT_SINGLE
#endif

/** Include mex source code for each template data type */

#define LTFAT_DOUBLE
#include "ltfat_mex_arg_helper.h"
#include MEX_FILE
#undef LTFAT_DOUBLE

#define LTFAT_SINGLE
#include "ltfat_mex_arg_helper.h"
#include MEX_FILE
#undef LTFAT_SINGLE

#include "mex.h"
#include <stdio.h>
#include <string.h>

#define LTFAT_MEXERRMSG(s,...)                 \
        char sChars[256];                      \
        snprintf(sChars,256,(s),__VA_ARGS__);  \
        mexErrMsgTxt(sChars)


/** C99 headers for a generic complex number manipulations
#include <complex.h>
#include <tgmath.h>
*/

/** Function prototypes */
bool checkIsSingle(const mxArray *prhsEl);
mxArray* recastToSingle(mxArray* prhsEl);
void checkArgs(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);

void checkArgs(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    #ifdef NARGINEQ
       if(nrhs!=NARGINEQ)
       {
          LTFAT_MEXERRMSG("Expected %i input arguments. Only %i passed.",NARGINEQ,nrhs);
       }
    #endif
    #if defined(NARGINLE)&&!defined(NARGINEQ)
       if(nrhs<=NARGINLE)
       {
          LTFAT_MEXERRMSG("Too many input arguments. Expected %i or less input arguments. Passed %i arg.",NARGINLE,nrhs);
       }
    #endif
    #if defined(NARGINGE)&&!defined(NARGINEQ)
       if(nrhs>=NARGINLE)
       {
          LTFAT_MEXERRMSG("Too few input arguments. Expected %i or more input arguments. Passed %i arg.",NARGINLE,nrhs);
       }
    #endif

}

/* Helper recasting functions */
mxArray* recastToSingle(mxArray* prhsEl)
{
   // return the input pointer if the input parameter alrady contains single prec. data
   if(checkIsSingle(prhsEl))
   {
     return prhsEl;
   }

   // if the input is cell array, cast all it's elements to single
   if(mxIsCell(prhsEl))
   {
      mxArray* tmpCell = mxCreateCellMatrix(mxGetM(prhsEl), mxGetN(prhsEl));
      for(unsigned int jj=0;jj<mxGetNumberOfElements(prhsEl);jj++)
      {
         mxSetCell(tmpCell, (mwIndex) jj, recastToSingle(mxGetCell(prhsEl, jj)));
      }
      return tmpCell;
   }

   if(mxIsStruct(prhsEl))
   {
      int nfields = mxGetNumberOfFields(prhsEl);
      const mwSize *dims = mxGetDimensions(prhsEl);
      mwSize ndim = mxGetNumberOfDimensions(prhsEl);

      const char **fieldnames = mxMalloc(nfields*sizeof(const char *));
      for(int ii=0;ii<nfields;ii++)
      {
         fieldnames[ii] = mxGetFieldNameByNumber(prhsEl,ii);
      }
      mxArray* tmpStructArr = mxCreateStructArray(ndim,dims,nfields,fieldnames);
      for(unsigned int jj=0;jj<mxGetNumberOfElements(prhsEl);jj++)
      {
         for(int ii=0;ii<nfields;ii++)
         {
            mxSetFieldByNumber(prhsEl,jj,ii,recastToSingle(mxGetFieldByNumber(prhsEl,jj,ii)));
         }
      }
         return tmpStructArr;
   }


   // just copy pointer if the element is not numeric
   if(!mxIsNumeric(prhsEl))
   {
      return prhsEl;
   }

   mwSize ndim = mxGetNumberOfDimensions(prhsEl);
   const mwSize* dims = mxGetDimensions(prhsEl);
   mxArray* tmpEl = 0;
   unsigned int elToCopy = mxGetNumberOfElements(prhsEl);
   if(mxIsComplex(prhsEl))
   {
      tmpEl = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxCOMPLEX);
      elToCopy *= 2;
   }
   else
   {
      tmpEl = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
   }

   double* prhsElPtr = (double*) mxGetData(prhsEl);
   float* tmpElPtr = (float*) mxGetData(tmpEl);
   for(unsigned int jj=0;jj<elToCopy;jj++)
   {
      *(tmpElPtr++) = (float)( *(prhsElPtr++) );
   }
   return tmpEl;
}

bool checkIsSingle(const mxArray *prhsEl)
{
   if(mxIsCell(prhsEl))
   {
      for(unsigned int jj=0;jj<mxGetNumberOfElements(prhsEl);jj++)
      {
         if(checkIsSingle(mxGetCell(prhsEl, jj)))
            return true;
      }
      return false;
   }

    if(mxIsSingle(prhsEl))
    {
       return true;
    }
    else
    {
       return false;
    }
}

/** MEX entry function
    Handles recasting all defined inputs to a defined data type
 */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
 {
    checkArgs(nlhs, plhs,nrhs,prhs);
    bool isAnySingle = false;
    #ifdef DATATYPECHECK
    const int PRHSTOCHECK[] = { DATATYPECHECK };
    unsigned int PRHSTOCHECKlen = sizeof(PRHSTOCHECK)/sizeof(PRHSTOCHECK[0]);
    for(unsigned int ii=0;ii<PRHSTOCHECKlen;ii++)
    {
      if(checkIsSingle(prhs[PRHSTOCHECK[ii]]))
      {
          isAnySingle = true;
          break;
      }
    }
    #endif


    if(isAnySingle)
    {
        const mxArray **prhsAlt = mxMalloc(nrhs*sizeof(mxArray *));
        memcpy((void *)prhsAlt,(void *)prhs,nrhs*sizeof(mxArray *));

        for(unsigned int ii =0;ii<PRHSTOCHECKlen;ii++)
        {
           prhsAlt[PRHSTOCHECK[ii]] = recastToSingle((mxArray *)prhs[PRHSTOCHECK[ii]]);
        }

        TEMPLATE(preMexFn,float)(nlhs, plhs, nrhs, prhsAlt);
        TEMPLATE(MEX_FUNC,float)(nlhs,plhs,nrhs, prhsAlt);
        TEMPLATE(postMexFn,float)(nlhs, plhs, nrhs, prhsAlt);
    }
    else
    {
        TEMPLATE(preMexFn,double)(nlhs, plhs, nrhs, prhs);
        TEMPLATE(MEX_FUNC,double)(nlhs,plhs,nrhs,prhs);
        TEMPLATE(postMexFn,double)(nlhs, plhs, nrhs, prhs);
    }
 }

#endif // _LTFAT_MEX_TEMPLATEHELPER_H
#endif

