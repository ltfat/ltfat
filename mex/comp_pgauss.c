#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_pgauss(L,w,c_t,c_f);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int L;
   double w, c_t, c_f;
   double *g;
   ltfat_complex *gc;
  
   L=(int)mxGetScalar(prhs[0]);
   w=(double)mxGetScalar(prhs[1]);
   c_t=(double)mxGetScalar(prhs[2]);
   c_f=(double)mxGetScalar(prhs[3]);

  if (c_f==0.0)
  {
     plhs[0] = mxCreateDoubleMatrix(L, 1, mxREAL);
     g = mxGetPr(plhs[0]);

     pgauss(L, w, c_t,(double*)g);
  }
  else
  {
    gc = mxMalloc(L*sizeof(ltfat_complex));

    pgauss_cmplx(L, w, c_t,c_f,gc);

    plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);

    combined2split(L, (const ltfat_complex*)gc, plhs[0]);

    mxFree(gc);
  }

  return;
  
}


