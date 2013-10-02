
#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

/** ISNARGINEQ, ISNARGINLE, ISNARGINGE
    AT COMPILE-TIME:
    AT RUNTIME: Ensures correct number of the input parameters.
    WHEN MISSING: No input argument checks ale included in the final code.
*/
#define ISNARGINEQ 5
/** TYPEDEPARGS
    AT COMPILE-TIME: Defines integer array from the specified values.
    AT RUNTIME: The array is used to identify input arguments to be checked/reformated. Accepted inputs are numeric arrays,
                cell arrays containing only numeric arrays, structures having at least one field beeing numeric array.
    WHEN MISSING: No input modifications/checks are included in the code.
*/
#define TYPEDEPARGS 0, 1
/** SINGLEARGS
    AT COMPILE-TIME: Includes this file for the second time with TYPEDEPARGS input args. recast to float arrays (cells, structs).
    AT RUNTIME: If at least one of the TYPEDEPARGS input args. is float (single in MatLab), all TYPEDEPARGS are recast to floats.
    WHEN MISSING: TYPEDEPARGS input args can be only double arrays.
*/
#define SINGLEARGS
/** COMPLEXARGS, REALARGS
    AT COMPILE-TIME: (COMPLEXARGS) adds code for on-the-fly conversion from the Matlab complex number format to the
                     complex.h (interleaved) complex data format.
                     (REALARGS) and (COMPLEXARGS) allows both real and complex inputs. Have to be handled here.
    AT RUNTIME: (COMPLEXARGS) TYPEDEPARGS input args are recast to complex format even in they are real.
                (REALARGS) TYPEDEPARGS args are accepted only if they are real.
                (REALARGS) and (COMPLEXARGS) If at least one of the TYPEDEPARGS is complex do as (COMPLEXARGS), otherwise let
                the inputs untouched.
    WHEN MISSING: Real/Complex are not checked. No complex data format change.
*/

/** COMPLEXINDEPENDENT
    AT COMPILE-TIME: As if both COMPLEXARGS, REALARGS were defined.
    AT RUNTIME: As if both COMPLEXARGS, REALARGS were defined plus it is assumed that the called functions from the LTFAT
                backend are from ltfat_typecomplexindependent.h, e.i. there are
    WHEN MISSING: No input checks REAL/COMPLEX checks are included in the final code.
*/
#define COMPLEXINDEPENDENT
#define EXPORTALIAS comp_filterbank_td

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

/* Obtain this filename. */
#if defined(__GNUC__) || defined(__ICC)
  #define MEX_FILE __BASE_FILE__
#endif
#include "ltfat_mex_template_helper.h"


#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
/** USER DEFINED HEADERS **/
#include "math.h"

/*
%COMP_FILTERBANK_TD   Non-uniform filterbank
%   Usage:  c=comp_filterbank_td(f,g,a,offset,ext);
%
%   Input parameters:
%         f   : Input data - L*W array.
%         g   : Filterbank filters - length M cell-array of vectors of lengths filtLen(m).
%         a   : Subsampling factors - array of length M.
%         offset: Offset of the filters - scalar or array of length M.
%         ext : Border exension technique.
%
%   Output parameters:
%         c  : Cell array of length M. Each element is N(m)*W array.
*/

void LTFAT_NAME(tdMexAtExitFnc)()
{
   #ifdef _DEBUG
   mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
   #endif
}


void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  #ifdef _DEBUG
  static int atExitFncRegistered = 0;
  if(!atExitFncRegistered)
  {
     LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(tdMexAtExitFnc));
     atExitFncRegistered = 1;
  }
  #endif

  const mxArray* mxf = prhs[0];
  const mxArray* mxg = prhs[1];
  double* a = (double*) mxGetData(prhs[2]);
  double* offset = (double*) mxGetData(prhs[3]);
  char* ext = mxArrayToString(prhs[4]);


  // input data length
  mwSize L = mxGetM(mxf);
  // number of channels
  mwSize W = mxGetN(mxf);
  // filter number
  mwSize M = mxGetNumberOfElements(mxg);
  // filter lengths
  mwSize filtLen[M];
  //mwSize* filtLen = mxMalloc(M*sizeof(mwSize));
  for(unsigned int m=0;m<M;m++)
  {
     filtLen[m] = (mwSize) mxGetNumberOfElements(mxGetCell(mxg,m));
  }

  // output lengths
  mwSize outLen[M];
  //mwSize* outLen = mxMalloc(M*sizeof(mwSize));
  if(!strcmp(ext,"per"))
  {
     for(unsigned int m = 0; m < M; m++)
     {
        outLen[m] = (mwSize) ceil( L/a[m] );
     }
  }
  else if(!strcmp(ext,"valid"))
  {
     for(unsigned int m = 0; m < M; m++)
     {
        outLen[m] = (mwSize) ceil( (L-(filtLen[m]-1))/a[m] );
     }
  }
  else
  {
     for(unsigned int m = 0; m < M; m++)
     {
        outLen[m] = (mwSize) ceil( (L + filtLen[m] - 1 + offset[m] )/a[m] );
     }
  }

     // POINTER TO THE INPUT
     LTFAT_TYPE* fPtr = (LTFAT_TYPE*) mxGetData(prhs[0]);

     // POINTER TO THE FILTERS
     LTFAT_TYPE* gPtrs[M];
     // LTFAT_TYPE** gPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     for(mwIndex m=0;m<M;m++)
     {
        gPtrs[m] = (LTFAT_TYPE*) mxGetPr(mxGetCell(mxg, m));
     }

     // POINTER TO OUTPUTS
     LTFAT_TYPE* cPtrs[M]; // C99 feature
     //LTFAT_TYPE** cPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     plhs[0] = mxCreateCellMatrix(M, 1);
     for(mwIndex m=0;m<M;++m)
     {
        mxSetCell(plhs[0], m, ltfatCreateMatrix(outLen[m], W,LTFAT_MX_CLASSID,LTFAT_MX_COMPLEXITY));
        cPtrs[m] = (LTFAT_TYPE*) mxGetData(mxGetCell(plhs[0],m));
        memset(cPtrs[m],0,outLen[m]*W*sizeof(LTFAT_TYPE));
     }

     // over all channels
   //  #pragma omp parallel for private(m)

        for(mwIndex m =0; m<M; m++)
        {
          for(mwIndex w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_TYPE *fPtrCol = fPtr + w*L;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_TYPE *cPtrCol = cPtrs[m] + w*outLen[m];
           //conv_td_sub(fPtrCol,L,&cPtrCol,outLen[m],(const double**)&gPtrs[m],filtLen[m],1,a[m],skip[m],ltfatExtStringToEnum(ext),0);
           LTFAT_NAME(convsub_td)(fPtrCol,L,cPtrCol,outLen[m],gPtrs[m],filtLen[m],a[m],-offset[m],ltfatExtStringToEnum(ext));
          }
        }
}
#endif

