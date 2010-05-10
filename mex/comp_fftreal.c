#include "mex.h"
#include "config.h"
#include "fftw3.h"

/* Calling convention:
 *  comp_fftreal(f);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 

  int L, W, L2;
  double *f, *cout_r, *cout_i;
  fftw_iodim dims[1], howmanydims[1];
  fftw_plan p;

  L = mxGetM(prhs[0]); 
  W = mxGetN(prhs[0]);

  L2 = (L/2)+1;

  /* Get pointer to input. */
  f=mxGetPr(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(L2, W, mxCOMPLEX);

  /* Get pointer to output. */
  cout_r = mxGetPr(plhs[0]);
  cout_i = mxGetPi(plhs[0]);  

  /* Create plan. Copy data from f to cout. */
  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L;
  howmanydims[0].os = L2;

  /* The calling prototype
  fftw_plan fftw_plan_guru_split_dft_r2c(
          int rank, const fftw_iodim *dims,
          int howmany_rank, const fftw_iodim *howmany_dims,
          double *in, double *ro, double *io,
          unsigned flags);
  */

  p = fftw_plan_guru_split_dft_r2c(1, dims,
				   1, howmanydims,
				   f, cout_r, cout_i, FFTW_OPTITYPE);
  
  /* Real FFT. */
  fftw_execute(p);   
  
  fftw_destroy_plan(p);   
  
  return;  

}
