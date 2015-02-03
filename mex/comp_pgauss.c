/*
#if defined(_WIN32) || defined(__WIN32__)
#  define DLL_EXPORT_SYM __declspec(dllexport)
#else
#  define EXPORT_EXTERN_C __attribute__((visibility("default")))
#endif

#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"
*/


#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_DOUBLE)
#include "ltfat_types.h"


/* Calling convention:
 *  comp_pgauss(L,w,c_t,c_f);
 */

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
		                        int UNUSED(nrhs), const mxArray *prhs[] )

{
   int L;
   double w, c_t, c_f;
   double *g;

   L=(int)mxGetScalar(prhs[0]);
   w=(double)mxGetScalar(prhs[1]);
   c_t=(double)mxGetScalar(prhs[2]);
   c_f=(double)mxGetScalar(prhs[3]);

  if (c_f==0.0)
  {
     plhs[0] = mxCreateDoubleMatrix(L, 1, mxREAL);
     g = mxGetPr(plhs[0]);

     pgauss_d(L, w, c_t,(double*)g);
  }
  else
  {
    plhs[0] = ltfatCreateMatrix(L, 1, mxDOUBLE_CLASS, mxCOMPLEX);
    LTFAT_COMPLEX *gc =  mxGetData(plhs[0]);
    pgauss_cmplx_d(L, w, c_t,c_f,gc);
  }

  return;

}

#endif


