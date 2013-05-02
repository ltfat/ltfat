//#define MEX_FILE
/** Allow including this file only if MEX_FILE is defined */
#if defined(MEX_FILE)

/** Allow including this file only once */
#ifndef _LTFAT_MEX_TEMPLATEHELPER_H
#define _LTFAT_MEX_TEMPLATEHELPER_H 1

/** Adds symbol exporting function decorator to mexFunction.
    On windows, def file is no longer needed. For MinGW, it
    suppresses the default "export-all-symbols" behavior. **/
#if defined(_WIN32) || defined(__WIN32__)
#  define DLL_EXPORT_SYM __declspec(dllexport)
#endif

/** Template macros */
#define LTFAT_CAT(prefix,name) prefix##name
#define LTFAT_TEMPLATE(prefix,name) LTFAT_CAT(prefix,name)

/** Undefine LTFAT_DOUBLE, LTFAT_SINGLE and LTFAT_COMPLEXTYPE if they are set */
#ifdef LTFAT_DOUBLE
#  undef LTFAT_DOUBLE
#endif
#ifdef LTFAT_SINGLE
#  undef LTFAT_SINGLE
#endif
#ifdef LTFAT_COMPLEXTYPE
#  undef LTFAT_COMPLEXTYPE
#endif



#include "ltfat.h"
#include <stdio.h>
#include <string.h>
#include "mex.h"
/** C99 headers for a generic complex number manipulations */
#if (defined(COMPLEXINDEPENDENT)||defined(COMPLEXARGS)) && !defined(NOCOMPLEXFMTCHANGE)
#  include <complex.h>
//#  include <tgmath.h>
#endif

/** Helper MACROS */
#ifdef _DEBUG
#define DEBUGINFO  mexPrintf("File: %s, func: %s \n",__BASE_FILE__,__func__);
#else
#define DEBUGINFO
#endif


/** Helper function headers, to allow them to be used in the MEX_FILE */
mxArray *ltfatCreateMatrix(mwSize M, mwSize N,mxClassID classid,mxComplexity complexFlag);
mxArray *ltfatCreateNdimArray(mwSize ndim, const mwSize *dims,mxClassID classid,mxComplexity complexFlag);

/** Include mex source code for each template data type */

#define LTFAT_DOUBLE
#include "ltfat_mex_typeindependent.h"
#include "ltfat_mex_typecomplexindependent.h"
#include MEX_FILE
#ifdef COMPLEXINDEPENDENT
#  define LTFAT_COMPLEXTYPE
#  include "ltfat_mex_typecomplexindependent.h"
#  include MEX_FILE
#  undef LTFAT_COMPLEXTYPE
#endif
#undef LTFAT_DOUBLE

#ifdef SINGLEARGS
#  define LTFAT_SINGLE
#  include "ltfat_mex_typeindependent.h"
#  include "ltfat_mex_typecomplexindependent.h"
#  include MEX_FILE
#  ifdef COMPLEXINDEPENDENT
#    define LTFAT_COMPLEXTYPE
#    include "ltfat_mex_typecomplexindependent.h"
#    include MEX_FILE
#    undef LTFAT_COMPLEXTYPE
#  endif
#  undef LTFAT_SINGLE
#endif

/**
* MACRO-FU
*/
#define LTFAT_MEXERRMSG(s,...)                 \
        do{                                    \
        char sChars[256];                      \
        snprintf(sChars,255,(s),__VA_ARGS__);  \
        mexErrMsgTxt(sChars);                  \
        }while(0)

#define FOREACH(item, array)                            \
        for(size_t keep = 1,                            \
            count = 0,                                  \
            size = sizeof (array) / sizeof *(array);    \
        keep && count != size;                          \
        keep = !keep, count++)                          \
        for(item = (array) + count; keep; keep = !keep)


#define FORSUBSET(item, array, subset)                            \
        for(size_t keep = 1,                                      \
            count = 0,                                            \
            size = sizeof (subset) / sizeof *(subset);            \
        keep && count != size;                                    \
        keep = !keep, count++)                                    \
        for(item = (array) + (subset)[count]; keep; keep = !keep)

#define FORSUBSETIDX(itemIdx, array, subset)                 \
        for(size_t keep = 1,                                      \
            count = 0,                                            \
            size = sizeof (subset) / sizeof *(subset);            \
        keep && count != size;                                    \
        keep = !keep, count++)                                    \
        for(itemIdx=&(subset)[count]; keep; keep = !keep)


/** Function prototypes */
bool checkIsReal(const mxArray *prhsEl);
bool checkIsSingle(const mxArray *prhsEl);
mxArray* recastToSingle(mxArray* prhsEl);
void checkArgs(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);
mwSize sizeofClassid(mxClassID classid);


mwSize sizeofClassid(mxClassID classid)
{
 switch(classid)
 {
    case mxSINGLE_CLASS:
       return sizeof(float);
    case mxDOUBLE_CLASS:
       return sizeof(double);
    case mxUNKNOWN_CLASS:
    case mxCELL_CLASS:
    case mxSTRUCT_CLASS:
    case mxFUNCTION_CLASS:
    default:
       mexErrMsgTxt("Usnupported data type. Add more if you need..");
       return 0;
    break;
 }
}

mxArray *ltfatCreateMatrix(mwSize M, mwSize N,mxClassID classid,mxComplexity complexFlag)
{
  const mwSize dims[]= {M,N};
  return ltfatCreateNdimArray(2, dims,classid,complexFlag);
}

mxArray *ltfatCreateNdimArray(mwSize ndim, const mwSize *dims,mxClassID classid, mxComplexity complexFlag)
{
   if(complexFlag==mxREAL)
     return mxCreateNumericArray(ndim,dims,classid,mxREAL);

   #if (!defined(COMPLEXINDEPENDENT) && !defined(COMPLEXARGS) && !defined(REALARGS)) || defined(NOCOMPLEXFMTCHANGE)
   if(complexFlag==mxCOMPLEX)
     return mxCreateNumericArray(ndim,dims,classid,mxCOMPLEX);
   #else
   if(complexFlag==mxCOMPLEX)
   {
      // Ugly...
      mwIndex dummyndim = 1;
      const mwSize dummyDims[] = {1};
      mxArray* out = mxCreateNumericArray(dummyndim,dummyDims,classid,mxCOMPLEX);
      // Set correct dimensions
      mxSetDimensions(out,dims,ndim);
      mwSize L = mxGetNumberOfElements(out);

      mxFree(mxGetPr(out));
      mxFree(mxGetPi(out));
      mxSetData(out,(void*)mxCalloc(L,2*sizeofClassid(classid)));
      mxSetImagData(out,mxCalloc(1,1));
      // To avoid automatic deallocation by the MEX memory manager

      return out;
   }
   #endif


   mexErrMsgTxt("Error in ltfatCreateNumericArray. Very strange...");
   return NULL;
}

void checkArgs(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    #ifdef ISNARGINEQ
       if(nrhs!=ISNARGINEQ)
       {
          LTFAT_MEXERRMSG("Expected %i input arguments. Only %i passed.",ISNARGINEQ,nrhs);
       }
    #endif
    #if defined(ISNARGINLE)&&!defined(ISNARGINEQ)
       if(nrhs<=ISNARGINLE)
       {
          LTFAT_MEXERRMSG("Too many input arguments. Expected %i or less input arguments. Passed %i arg.",ISNARGINLE,nrhs);
       }
    #endif
    #if defined(ISNARGINGE)&&!defined(ISNARGINEQ)
       if(nrhs>=ISNARGINLE)
       {
          LTFAT_MEXERRMSG("Too few input arguments. Expected %i or more input arguments. Passed %i arg.",ISNARGINLE,nrhs);
       }
    #endif

}

/* Helper recasting functions */
mxArray* recastToSingle(mxArray* prhsEl)
{
   // return the input pointer if the input parameter already contains single prec. data
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
      mwSize nfields = mxGetNumberOfFields(prhsEl);
      const mwSize *dims = mxGetDimensions(prhsEl);
      mwSize ndim = mxGetNumberOfDimensions(prhsEl);

      const char **fieldnames = mxMalloc(nfields*sizeof(const char *));
      for(mwSize ii=0;ii<nfields;ii++)
      {
         fieldnames[ii] = mxGetFieldNameByNumber(prhsEl,ii);
      }
      mxArray* tmpStructArr = mxCreateStructArray(ndim,dims,nfields,fieldnames);
      for(mwIndex jj=0;jj<mxGetNumberOfElements(prhsEl);jj++)
      {
         for(mwIndex ii=0;ii<nfields;ii++)
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
   mwSize elToCopy = mxGetNumberOfElements(prhsEl);

   if(mxIsComplex(prhsEl))
   {
      tmpEl = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxCOMPLEX);
   }
   else
   {
      tmpEl = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
   }

   double* prhsElPtr = (double*) mxGetPr(prhsEl);
   float* tmpElPtr = (float*) mxGetPr(tmpEl);

   for(mwIndex jj=0;jj<elToCopy;jj++)
   {
      *(tmpElPtr++) = (float)( *(prhsElPtr++) );
   }

   if(mxIsComplex(prhsEl))
   {
      double* prhsElPtr_i = (double*) mxGetPi(prhsEl);
      float* tmpElPtr_i = (float*) mxGetPi(tmpEl);

      for(mwIndex jj=0;jj<elToCopy;jj++)
      {
         *(tmpElPtr_i++) = (float)( *(prhsElPtr_i++) );
      }
   }

   return tmpEl;
}

bool checkIsSingle(const mxArray *prhsEl)
{
   if(mxIsCell(prhsEl))
   {
      for(mwIndex jj=0;jj<mxGetNumberOfElements(prhsEl);jj++)
      {
         if(checkIsSingle(mxGetCell(prhsEl, jj)))
            return true;
      }
      return false;
   }
   return mxIsSingle(prhsEl);
}

bool checkIsReal(const mxArray *prhsEl)
{
   if(mxIsCell(prhsEl))
   {
      bool isAllReal = false;
      for(mwIndex jj=0;jj<mxGetNumberOfElements(prhsEl);jj++)
      {
         if(!(isAllReal = checkIsReal(mxGetCell(prhsEl, jj))))
            break;
      }
      return isAllReal;
   }

    return !mxIsComplex(prhsEl);
}

/** MEX entry function
    Handles recasting all defined inputs to a defined data type
 */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
 {
   #if defined(ISNARGINEQ) || defined(ISNARGINLE) || defined(ISNARGINGE)
     checkArgs(nlhs, plhs,nrhs,prhs);
   #endif
   #ifndef TYPEDEPARGS
       LTFAT_NAME_DOUBLE(ltfatMexFnc)(nlhs,plhs,nrhs,prhs); // pass trough
   #else
      // Create array of indexes to be checked
      int prhsToCheck[] = { TYPEDEPARGS };

      // Check if any of the indexes is less than nrhs
      int maxPrhsToCheck = 0;
      FOREACH(int* valPr,prhsToCheck)
         if(*valPr+1>maxPrhsToCheck)
            maxPrhsToCheck = *valPr;

      if(!((maxPrhsToCheck)<nrhs))
        LTFAT_MEXERRMSG("To few input arguments. Expected at least %i args.",maxPrhsToCheck+1);

      // Indicator array defining which input arg. should be reformated.
      bool recastToComplexIndArr[nrhs];
      memset(recastToComplexIndArr,0,sizeof(recastToComplexIndArr));

      bool isAnyComplex = false;
      FORSUBSET(const mxArray **prhsElPtr, prhs, prhsToCheck)
        if( (isAnyComplex = !checkIsReal(*prhsElPtr))) break;

      #if !defined(COMPLEXARGS) && defined(REALARGS) && !defined(COMPLEXINDEPENDENT)
      FORSUBSET(const mxArray **prhsElPtr, prhs, prhsToCheck)
         if( !checkIsReal(*prhsElPtr))
         {
            mexErrMsgTxt("Complex input arguments are not alowed.");
            break;
         }
      #endif
      #if (defined(COMPLEXARGS) && !defined(REALARGS)) && !defined(COMPLEXINDEPENDENT)
      FORSUBSETIDX(int*  prhsElIdx,prhs, prhsToCheck)
         recastToComplexIndArr[*prhsElIdx]=true;
      #endif

      #if (defined(COMPLEXARGS) && defined(REALARGS)) || defined(COMPLEXINDEPENDENT)
      if(isAnyComplex)
      FORSUBSETIDX(int*  prhsElIdx,prhs, prhsToCheck)
          recastToComplexIndArr[*prhsElIdx]=true;
      #endif

      const mxArray **prhsAlt = mxMalloc(nrhs*sizeof(mxArray *));
      memcpy((void *)prhsAlt,(void *)prhs,nrhs*sizeof(mxArray *));

      bool isAnySingle = false;
      #ifdef SINGLEARGS
      FORSUBSET(const mxArray **prhsElPtr, prhs, prhsToCheck)
        if( (isAnySingle = checkIsSingle(*prhsElPtr))) break;

      if(isAnySingle)
      {
        FORSUBSETIDX(int*  prhsElIdx,prhs, prhsToCheck)
           prhsAlt[*prhsElIdx] = recastToSingle((mxArray *)prhs[*prhsElIdx]);

        #if (defined(COMPLEXINDEPENDENT) || defined(COMPLEXARGS)) && !defined(NOCOMPLEXFMTCHANGE)
        for(int ii=0;ii<nrhs;ii++)
           if(recastToComplexIndArr[ii])
              prhsAlt[ii] = LTFAT_NAME_SINGLE(mexSplit2combined)(prhsAlt[ii]);
        #endif

        #if defined(COMPLEXINDEPENDENT)

        bool isAllReal = false;
        bool isAllComplex = false;

        FORSUBSET(const mxArray **prhsElPtr, prhsAlt, prhsToCheck)
          if( !(isAllReal = checkIsReal(*prhsElPtr))) break;

        FORSUBSET(const mxArray **prhsElPtr, prhsAlt, prhsToCheck)
          if( !(isAllComplex = !checkIsReal(*prhsElPtr))) break;

        if(!(isAllReal ^ isAllComplex))
           mexErrMsgTxt("Template subsystem error. My bad .");

        if(isAllReal)
           LTFAT_NAME_SINGLE(ltfatMexFnc)(nlhs,plhs,nrhs, prhsAlt);
        else if(isAllComplex)
           LTFAT_NAME_COMPLEXSINGLE(ltfatMexFnc)(nlhs,plhs,nrhs, prhsAlt);
        #else
        LTFAT_NAME_SINGLE(ltfatMexFnc)(nlhs,plhs,nrhs, prhsAlt);
        #endif //COMPLEXINDEPENDENT

        #if (defined(COMPLEXINDEPENDENT) || defined(COMPLEXARGS) || defined(REALARGS)) && !defined(NOCOMPLEXFMTCHANGE)
        if(!checkIsReal(plhs[0]))
        plhs[0] = LTFAT_NAME_SINGLE(mexCombined2split)(plhs[0]);
        #endif
      }
      #endif // SINGLEARGS
      if(!isAnySingle)
      {

        #if (defined(COMPLEXINDEPENDENT) || defined(COMPLEXARGS)) && !defined(NOCOMPLEXFMTCHANGE)
        for(int ii=0;ii<nrhs;ii++)
           if(recastToComplexIndArr[ii])
              prhsAlt[ii] = LTFAT_NAME_DOUBLE(mexSplit2combined)(prhsAlt[ii]);

        #endif


        #if defined(COMPLEXINDEPENDENT)

        bool isAllReal = false;
        bool isAllComplex = false;

        FORSUBSET(const mxArray **prhsElPtr, prhsAlt, prhsToCheck)
          if( !(isAllReal = checkIsReal(*prhsElPtr))) break;

        FORSUBSET(const mxArray **prhsElPtr, prhsAlt, prhsToCheck)
          if( !(isAllComplex = !checkIsReal(*prhsElPtr))) break;

        if((isAllReal == isAllComplex))
           mexErrMsgTxt("Template subsystem error. My bad...");

        if(isAllReal)
           LTFAT_NAME_DOUBLE(ltfatMexFnc)(nlhs,plhs,nrhs,prhsAlt);
        else if(isAllComplex)
           LTFAT_NAME_COMPLEXDOUBLE(ltfatMexFnc)(nlhs,plhs,nrhs, prhsAlt);

        #else
          LTFAT_NAME_DOUBLE(ltfatMexFnc)(nlhs,plhs,nrhs,prhsAlt);
        #endif


        #if (defined(COMPLEXINDEPENDENT) || defined(COMPLEXARGS) || defined(REALARGS)) && !defined(NOCOMPLEXFMTCHANGE)
        if(!checkIsReal(plhs[0]))
          plhs[0] = LTFAT_NAME_DOUBLE(mexCombined2split)(plhs[0]);

        #endif


      }



    #endif // TYPEDEPARGS
 }

#endif // _LTFAT_MEX_TEMPLATEHELPER_H
#endif // defined(MEX_FILE)


