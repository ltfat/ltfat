#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *    c=comp_nonsepdgt_multi(f,g,a,M,lt);
 */

int ltfat_round(double x)
{
  return (int)(x+.5); 
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int a, M, N, L, W, Lg, lt1, lt2;
   mwSize ndim;
   mwSize dims[3];

   ltfat_complex *f_combined, *g_combined, *out_combined;

   double *lt;


   /* Get matrix dimensions.*/
   L  = mxGetM(prhs[0]); 
   W  = mxGetN(prhs[0]); 
   Lg = mxGetM(prhs[1]); 
   
   a=(int)mxGetScalar(prhs[2]);
   M=(int)mxGetScalar(prhs[3]);

   /* Read the values of lt and round them to integers. */
   lt = mxGetPr(prhs[4]);
   lt1 = ltfat_round(lt[0]);
   lt2 = ltfat_round(lt[1]);

   N  = L/a;
  
   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   f_combined=(ltfat_complex*)mxMalloc(L*W*sizeof(ltfat_complex));
   g_combined=mxMalloc(Lg*sizeof(ltfat_complex));
   out_combined=(ltfat_complex*)mxMalloc(M*N*W*sizeof(ltfat_complex));

   split2combined(L*W, prhs[0], f_combined);
   split2combined(Lg, prhs[1], g_combined);

   dgt_multi((const ltfat_complex*)f_combined,
	     (const ltfat_complex*)g_combined,
	     L,Lg,W,a,M,lt1,lt2,
	     (ltfat_complex*)out_combined);
        
   mxFree(f_combined);   
   mxFree(g_combined);

   plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
      
   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);
   return;
}
