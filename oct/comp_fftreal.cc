#include <octave/oct.h>
#include "config.h"
#include "fftw3.h"

DEFUN_DLD (comp_fftreal, args, ,
  "This function calls the FFTW3 real FFT\n\
  c=comp_fftreal(f);\n\
  Yeah.")
{

  fftw_plan p;

  const Matrix f = args(0).matrix_value();

  const int L = f.rows();
  const int W = f.columns();
  
  const int L2 = (L/2)+1;
  
  ComplexMatrix cout(L2,W);  
  
  /* Create plan. Copy data from f to cout. */
  p = fftw_plan_many_dft_r2c(1, &L, W,
			     (double*)f.data(), NULL,
			     1, L,
			     (fftw_complex*)cout.data(), NULL,
			     1, L2,
			     FFTW_OPTITYPE);
  
  /* Real FFT. */
  fftw_execute(p);   
  
  fftw_destroy_plan(p);   
    
  return octave_value (cout);

}
