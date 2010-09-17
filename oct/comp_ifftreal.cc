#include <octave/oct.h>
#include "config.h"
#include "fftw3.h"

DEFUN_DLD (comp_ifftreal, args, ,
  "This function calls the FFTW3 real FFT\n\
  f=comp_fftreal(c,L2);\n\
  Yeah.")
{

  fftw_plan p;

  const ComplexMatrix c = args(0).complex_matrix_value();

  const int L2 = c.rows();
  const int W  = c.columns();
  
  const int L = args(1).int_value();
  
  Matrix f(L,W);  
  
  /* Prototype:
     fftw_plan fftw_plan_many_dft_c2r(int rank, const int *n, int howmany,
                                      fftw_complex *in, const int *inembed,
                                      int istride, int idist,
                                      double *out, const int *onembed,
                                      int ostride, int odist,
                                      unsigned flags);
  */

  /* Create plan. Copy data from f to cout. */
  p = fftw_plan_many_dft_c2r(1, &L, W,
			     (fftw_complex*)c.data(), NULL,
			     1, L2,
			     f.fortran_vec(), NULL,
			     1, L,
			     FFTW_ESTIMATE);
  
  /* Real IFFT. */
  fftw_execute(p);   

  /* Scale, because Octave's normalization is different. */
  double s  = 1.0/L;
  double *fp = f.fortran_vec();
  for (int ii=0; ii<L*W; ii++)
  {
      fp[ii] *=s;
  }
  fftw_destroy_plan(p);   
    
  return octave_value (f);

}
