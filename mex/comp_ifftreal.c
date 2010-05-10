#include "mex.h"
#include "config.h"
#include "fftw3.h"
#include "ltfat-mex-helper.h"

/* Calling convention:
 *  comp_ifftreal(f,N);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{   

  int ii, L, W, L2;
  fftw_plan p;
  double *f, s;

#ifdef FIXME_THIS_SECTION_CONTAINS_AN_ERROR 

  L2 = mxGetM(prhs[0]); 
  W  = mxGetN(prhs[0]);
  L  = (int)mxGetScalar(prhs[1]);

  /* Create output and get pointer */
  plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);
  f=mxGetPr(plhs[0]);

  /* This section is not being compiled. It contains a segmentation
   * faults. The idea is to pass Matlab's split memory layout directly
   * to FFTW
   */

  fftw_iodim dims[1], howmanydims[1];

  /* Create plan. Copy data from cin to f. */
  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L2;
  howmanydims[0].os = L;

  /* The calling prototype

    fftw_plan fftw_plan_guru_split_dft_c2r(
          int rank, const fftw_iodim *dims,
          int howmany_rank, const fftw_iodim *howmany_dims,
          double *ri, double *ii, double *out,
          unsigned flags);
  */

  p = fftw_plan_guru_split_dft_r2c(1, dims,
				   1, howmanydims,
				   mxGetPr(prhs[0]), mxGetPi(prhs[0]), 
				   f, FFTW_ESTIMATE);
  
  /* Real IFFT. */
  fftw_execute(p);   
    
#else

  /* This section copies to a combined layout, and use the exact same
   * FFTW call as the Octave interface. */

  fftw_complex *c_combined;

  L2 = mxGetM(prhs[0]); 
  W  = mxGetN(prhs[0]);
  L  = (int)mxGetScalar(prhs[1]);

  /* Create output and get pointer */
  plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);
  f=mxGetPr(plhs[0]);
  
  c_combined = mxMalloc(L2*W*sizeof(fftw_complex));

  split2combined(L2*W, prhs[0], c_combined);  

  /* Create plan. Copy data from f to cout. */
  p = fftw_plan_many_dft_c2r(1, &L, W,
			     c_combined, NULL,
			     1, L2,
			     f, NULL,
			     1, L,
			     FFTW_ESTIMATE);
  
  /* Real IFFT. */
  fftw_execute(p);   

  mxFree(c_combined);
  
#endif

  fftw_destroy_plan(p);     

  /* Scale, because FFTW's normalization is different. */
  s  = 1.0/L;
  for (ii=0; ii<L*W; ii++)
    {
      f[ii] *=s;
  }


  return;  

}
