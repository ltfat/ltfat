#include "mex.h"
#include "config.h"
#include "ltfat.h"

/* Calling convention:
 *  comp_dgt_fb(f,g,a,M);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
   
{ 
   int L, gl,W, a, M, N;
   double *f_combined, *g_combined, *out_combined;
   double *f_r,*f_i,*g_r,*g_i,*out_r,*out_i;
   
   int ii;
   
   /* Get matrix dimensions.*/
   L  = mxGetM(prhs[0]); 
   W  = mxGetN(prhs[0]);
   gl = mxGetM(prhs[1]); 

   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);
   
   N=L/a;
   
   /* Create temporary matrices to convert to correct complex layout. */
   
   f_combined=mxCalloc(L*W*2,sizeof(double));
   g_combined=mxCalloc(L*2,sizeof(double));
   out_combined=mxCalloc(M*N*W*2,sizeof(double));
   
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
   
   
   g_r=mxGetPr(prhs[1]);
   for (ii=0;ii<gl; ii++)
   {
      g_combined[2*ii]=g_r[ii];
   }

   if (mxIsComplex(prhs[1]))
   {
      g_i=mxGetPi(prhs[1]);
      for (ii=0;ii<gl; ii++)
      {
	 g_combined[2*ii+1]=g_i[ii];
      }
   }
   
   dgt_fb((const ltfat_complex*)f_combined,(const ltfat_complex*)g_combined,L,gl,W,a,M,
              (ltfat_complex*)out_combined);
   
   mxFree(f_combined);
   mxFree(g_combined);
   
   plhs[0] = mxCreateDoubleMatrix(M, N*W, mxCOMPLEX);
   
   out_r=mxGetPr(plhs[0]);
   out_i=mxGetPi(plhs[0]);
   
   for (ii=0;ii<M*N*W; ii++)
   {
      out_r[ii]=out_combined[2*ii];
      out_i[ii]=out_combined[2*ii+1];
   }
   
   return;
   
}


