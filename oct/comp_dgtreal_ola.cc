#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_dgtreal_ola, args, ,
  "This function calls the C-library\n\
  c=comp_dgtreal_ola(f,g,a,M,bl);\n\
  Yeah.")
{

  const int a  = args(2).int_value();
  const int M  = args(3).int_value();
  const int bl = args(4).int_value();

  const Matrix f = args(0).matrix_value();
  const Matrix g = args(1).matrix_value();
  
  const int L  = f.rows();
  const int W  = f.columns();
  const int gl = g.rows();
  const int N = L/a;

  const int M2 = M/2+1;
  
  dim_vector dims_out(M2,N,W);  
  dims_out.chop_trailing_singletons();
  
  ComplexNDArray cout(dims_out);  
  
  dgtreal_ola(f.data(),g.data(),
	  L,gl,W,a,M,bl,
	  (ltfat_complex*)cout.data());  
  
  return octave_value (cout);

}
