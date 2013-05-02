#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 1
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
/*
  Defining forwarders since there is no simple way of unified call to
  col2diag since the real and imag planes are processed separatelly
*/
#if defined(LTFAT_DOUBLE)
static inline void LTFAT_NAME(fwd_col2diag)(const double* cin, const int L, double* cout)
{
   col2diag_d(cin,L,cout);
}
#elif defined(LTFAT_SINGLE)
static inline void LTFAT_NAME(fwd_col2diag)(const float* cin, const int L, float* cout)
{
   col2diag_s(cin,L,cout);
}
#endif

/* Calling convention:
 *  cout=comp_col2diag(cin);
 */
void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    int L = mxGetM(prhs[0]);
    plhs[0] = ltfatCreateMatrix(L, L, LTFAT_MX_CLASSID, LTFAT_MX_COMPLEXITY);

    #ifdef NOCOMPLEXFMTCHANGE
    LTFAT_REAL* cout_r = (LTFAT_REAL*) mxGetPr(plhs[0]);
    LTFAT_REAL* cin_r = (LTFAT_REAL*) mxGetPr(prhs[0]);
    LTFAT_NAME(fwd_col2diag)(cin_r,L,cout_r);

      #ifdef LTFAT_COMPLEXTYPE
      // Treat the imaginary part
      LTFAT_REAL* cin_i= (LTFAT_REAL*) mxGetPi(prhs[0]);
      LTFAT_REAL* cout_i= (LTFAT_REAL*) mxGetPi(plhs[0]);
      LTFAT_NAME(fwd_col2diag)(cin_i, L,cout_i);
      #endif
    #else
      LTFAT_TYPE* cin_r= (LTFAT_TYPE*) mxGetData(prhs[0]);
      LTFAT_TYPE* cout_r= (LTFAT_TYPE*) mxGetData(plhs[0]);
      LTFAT_NAME(col2diag)(cin_r, L,cout_r);
    #endif
}
#endif




/* Calling convention:
 *  cout=comp_col2diag(cin);
 */
/*
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{


   int L;
   double *cin_r,*cin_i, *cout_r,*cout_i;

   // Get matrix dimensions.
   L = mxGetM(prhs[0]);

   cin_r=mxGetPr(prhs[0]);

   if (mxIsComplex(prhs[0]))
   {
      plhs[0] = mxCreateDoubleMatrix(L,L, mxCOMPLEX);
   }
   else
   {
      plhs[0] = mxCreateDoubleMatrix(L,L, mxREAL);
   }

   // Treat the real part
   cout_r=mxGetPr(plhs[0]);
   col2diag_r((double*)cin_r, L, (double*)cout_r);

   if (mxIsComplex(prhs[0]))
   {
      // Treat the imaginary part
      cin_i=mxGetPi(prhs[0]);
      cout_i=mxGetPi(plhs[0]);
      col2diag_r((double*)cin_i, L, (double*)cout_i);
   }

   return;

}
*/


