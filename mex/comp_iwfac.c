#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_iwfac(gf,L,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 

  int L, a, M, N;
  ltfat_complex *g_combined, *gf_combined;

  L=(int)mxGetScalar(prhs[1]);
  a=(int)mxGetScalar(prhs[2]);
  M=(int)mxGetScalar(prhs[3]);
  
  /* Create temporary matrices to convert to correct complex layout. */

  g_combined=mxMalloc(L*sizeof(ltfat_complex));
  gf_combined=mxMalloc(L*sizeof(ltfat_complex));

  /* Copy the data. */

  split2combined(L, prhs[0], gf_combined);

  iwfac((const ltfat_complex*)gf_combined, L, a, M, g_combined);

  plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);

  combined2split(L, (const ltfat_complex*)g_combined, plhs[0]);

  mxFree(gf_combined);
  mxFree(g_combined);

  return;
  
}


