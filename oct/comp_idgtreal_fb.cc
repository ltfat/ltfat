#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_idgtreal_fb, args, ,
  "This function calls the C-library\n\
  c=comp_idgt_fb(c,g,L,a,M);\n\
  Yeah.\n")
{

  const ComplexMatrix coef = args(0).complex_matrix_value();
  const Matrix g = args(1).matrix_value();
  const int L = args(2).int_value();
  const int a = args(3).int_value();
  const int M = args(4).int_value();
  const int N = L/a;
  const int W = coef.columns()/N;
  const int gl = g.rows();

  Matrix f(L,W);  
  
  idgtreal_fb((ltfat_complex*)coef.data(),g.data(),L,gl,W,a,M,
	  f.fortran_vec());

  return octave_value (f);
}
