#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2

#define TYPEDEPARGS 0
#define SINGLEARGS
#define REALARGS
#define COMPLEXARGS
#define NOCOMPLEXFMTCHANGE

#endif /* _LTFAT_MEX_FILE */

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "config.h"


static LTFAT_FFTW(plan) LTFAT_NAME(p_old) = 0;

static void LTFAT_NAME(dctMexAtExitFnc)()
{
   if (LTFAT_NAME(p_old) != 0)
   {
      LTFAT_FFTW(destroy_plan)(LTFAT_NAME(p_old));
   }
}


void
LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                         int UNUSED(nrhs), const mxArray *prhs[] )
{
   // Register exit function only once
   static int atExitFncRegistered = 0;
   if (!atExitFncRegistered)
   {
      LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(dctMexAtExitFnc));
      atExitFncRegistered = 1;
   }

   LTFAT_REAL *c_r, *c_i = NULL;
   const LTFAT_REAL *f_r, *f_i = NULL;
   dct_kind kind = DSTI;

   mwIndex L = mxGetM(prhs[0]);
   mwIndex W = mxGetN(prhs[0]);
   mwIndex type = (mwIndex) mxGetScalar(prhs[1]);

// Copy inputs and get pointers
   if ( mxIsComplex(prhs[0]))
   {
      f_i =  mxGetImagData(prhs[0]);
      plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID, mxCOMPLEX);
      c_i = mxGetImagData(plhs[0]);
   }
   else
   {
      plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID, mxREAL);
   }

   f_r = mxGetData(prhs[0]);
   c_r = mxGetData(plhs[0]);

   switch (type)
   {
   case 1:
      kind = DSTI;
      break;
   case 2:
      kind = DSTII;
      break;
   case 3:
      kind = DSTIII;
      break;
   case 4:
      kind = DSTIV;
      break;
   default:
      mexErrMsgTxt("Unknown type.");
   }



   LTFAT_FFTW(plan) p = LTFAT_NAME(dst_init)( L, W, c_r, kind);



   LTFAT_NAME(dctMexAtExitFnc)();
   LTFAT_NAME(p_old) = p;


   LTFAT_NAME(dst_execute)(p, f_r, L, W, c_r, kind);
   if ( mxIsComplex(prhs[0]))
   {
      LTFAT_NAME(dst_execute)(p, f_i, L, W, c_i, kind);
   }


   return;
}
#endif

