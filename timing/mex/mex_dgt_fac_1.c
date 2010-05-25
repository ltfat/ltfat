#include "mex.h"
#include "config.h"
#include "../dgt_fac_1.c"
#include "../gcd.c"

/* Calling convention:
 *  comp_dgt_fac(f,gf,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, W, R, a, M, N;
   double *f_combined, *gf_combined, *out_combined;
   double *f_r,*f_i,*gf_r,*gf_i,*out_r,*out_i;
   
   int ii;
   
   /* Get matrix dimensions.*/
   L = mxGetM(prhs[0]); 
   W = mxGetN(prhs[0]);
   R = mxGetM(prhs[1])*mxGetN(prhs[1])/L; 
 
   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   
   N=L/a;
   
   /* Create temporary matrices to convert to correct complex layout. */
   
   f_combined=mxCalloc(L*W*2,sizeof(double));
   gf_combined=mxCalloc(L*R*2,sizeof(double));
   out_combined=mxCalloc(M*N*R*W*2,sizeof(double));
   
   /* Copy the data. */
   
   f_r=mxGetPr(prhs[0]);
   for (ii=0;ii<L*W; ii++)
   {
      f_combined[2*ii]=f_r[ii];
   }
   
   if (mxIsComplex(prhs[0]))
   {
      f_i=mxGetPi(prhs[0]);
      for (ii=0;ii<L*W; ii++)
      {
	 f_combined[2*ii+1]=f_i[ii];
      }
   }
   
   
   gf_r=mxGetPr(prhs[1]);
   for (ii=0;ii<L*R; ii++)
   {
      gf_combined[2*ii]=gf_r[ii];
   }

   if (mxIsComplex(prhs[1]))
   {
      gf_i=mxGetPi(prhs[1]);
      for (ii=0;ii<L*R; ii++)
      {
	 gf_combined[2*ii+1]=gf_i[ii];
      }
   }
   
   dgt_fac_1((ltfat_complex*)f_combined,(ltfat_complex*)gf_combined,L,W,R,a,M,
	   (ltfat_complex*)out_combined);
   
   mxFree(f_combined);
   mxFree(gf_combined);
   
   plhs[0] = mxCreateDoubleMatrix(M, N*W*R, mxCOMPLEX);
   
   out_r=mxGetPr(plhs[0]);
   out_i=mxGetPi(plhs[0]);
   
   for (ii=0;ii<M*N*W*R; ii++)
   {
      out_r[ii]=out_combined[2*ii];
      out_i[ii]=out_combined[2*ii+1];
   }
   
   return;
   
}


