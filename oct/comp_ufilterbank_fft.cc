#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_ufilterbank_fft, args, ,
  "This function calls the C-library\n\
  c=comp_pfilt_fir(f,g,a);\n\
  Yeah.")
{

  const ComplexMatrix f = args(0).complex_matrix_value();
  const ComplexMatrix g = args(1).complex_matrix_value();
  const int a = args(2).int_value();
     
  const int L  = f.rows();
  const int W  = f.columns();
  const int gl = g.rows();
  const int M  = g.columns();
  const int N = L/a;
  
  dim_vector dims_out(N,M,W);  
  dims_out.chop_trailing_singletons();
  
  ComplexNDArray cout(dims_out);  

  ufilterbank_fft((const ltfat_complex*)f.data(),(const ltfat_complex*)g.data(),
		  L,gl,W,a,M,(ltfat_complex*)cout.data());  
  
  return octave_value (cout);

}
