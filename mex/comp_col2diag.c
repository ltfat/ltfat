#include "mex.h"
#include "config.h"
#include "ltfat.h"

/* Calling convention:
 *  cout=comp_col2diag(cin);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
   int L;
   double *cin_r,*cin_i, *cout_r,*cout_i;
   
   int ii;
   
   /* Get matrix dimensions.*/
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

   /* Treat the real part */
   cout_r=mxGetPr(plhs[0]);
   col2diag_r((double*)cin_r, L, (double*)cout_r);

   if (mxIsComplex(prhs[0]))
   {
      /* Treat the imaginary part */
      cin_i=mxGetPi(prhs[0]);   
      cout_i=mxGetPi(plhs[0]);
      col2diag_r((double*)cin_i, L, (double*)cout_i);
   }
               
   return;
  
}


