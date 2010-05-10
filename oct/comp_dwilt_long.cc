#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"


DEFUN_DLD (comp_dwilt_long, args, ,
  "This function calls the C-library\n\
  coef=comp_dwilt_long(f,g,M,L);\n\
  Yeah.")
{
  const int M = args(2).int_value();
  const int L = args(3).int_value();
  const int N = L/M;

  if (args(0).is_complex_type())
  {

    const ComplexMatrix f = args(0).complex_matrix_value();
    const ComplexMatrix g = args(1).complex_matrix_value();
    
    const int W = f.columns();

    dim_vector dims_out(2*M,N/2,W);  
    dims_out.chop_trailing_singletons();
    
    ComplexNDArray cout(dims_out);  
    
    dwilt_long((ltfat_complex*)f.data(),(ltfat_complex*)g.data(), L, W, M,
	       (ltfat_complex*)cout.data());

    return octave_value (cout);

  }
  else
  {

    const Matrix f = args(0).matrix_value();
    const Matrix g = args(1).matrix_value();
    
    const int W = f.columns();

    dim_vector dims_out(2*M,N/2,W);  
    dims_out.chop_trailing_singletons();
    
    NDArray cout(dims_out);  
    
    dwiltreal_long((double*)f.data(),(double*)g.data(), L, W, M,
		  (double*)cout.data());
	
    return octave_value (cout);

  }

}
