#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_dgt_ola, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_ola(f,g,a,M,bl);\n\
  Yeah.")
{

  const int a  = args(2).int_value();
  const int M  = args(3).int_value();
  const int bl = args(4).int_value();

  const ComplexMatrix f = args(0).complex_matrix_value();
  const ComplexMatrix g = args(1).complex_matrix_value();
  
  const int L  = f.rows();
  const int W  = f.columns();
  const int gl = g.rows();
  const int N = L/a;
  
  dim_vector dims_out(M,N,W);  
  dims_out.chop_trailing_singletons();
  
  ComplexNDArray cout(dims_out);  
  
  dgt_ola((ltfat_complex*)f.data(),(ltfat_complex*)g.data(),
	  L,gl,W,a,M,bl,
	  (ltfat_complex*)cout.data());  
  
  return octave_value (cout);

}
