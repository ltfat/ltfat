#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_wfac(g,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int L, R, a, M, N, c, d, p, q,h_a,h_m;
   ltfat_complex *g_combined, *gf_combined;
   
   /* Get matrix dimensions.*/
   L = mxGetM(prhs[0]); 
   R = mxGetN(prhs[0]); 
   
   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);
   
   N=L/a;
   
   c=gcd(a, M,&h_a, &h_m);
   p=a/c;
   q=M/c;
   d=N/q;

   /* Create temporary matrices to convert to correct complex layout. */   
   gf_combined=mxMalloc(L*R*sizeof(ltfat_complex));
   
   if (mxIsComplex(prhs[0]))
   {
      /* g is complex */
      g_combined=mxMalloc(L*R*sizeof(ltfat_complex));

      split2combined(L*R, prhs[0], g_combined);

      wfac((const ltfat_complex*)g_combined, 
	   L, R, a, M, (ltfat_complex*)gf_combined);

      mxFree(g_combined);
   }
   else
   {
      /* g is real */
     wfac_r(mxGetPr(prhs[0]), L, R, a, M, 
	    gf_combined);      
   }
               
   plhs[0] = mxCreateDoubleMatrix(p*q*R, c*d, mxCOMPLEX);
   
   combined2split(L*R, (const ltfat_complex*)gf_combined, plhs[0]);

   mxFree(gf_combined);

   return;
  
}


