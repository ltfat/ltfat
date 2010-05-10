#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_idgtreal_fac, args, ,
  "This function calls the C-library\n\
  c=comp_idgt_fac(c,gf,L,a,M);\n\
  Yeah.\n")
{

  const ComplexMatrix coef = args(0).complex_matrix_value();
  const ComplexMatrix gf = args(1).complex_matrix_value();
  const int L = args(2).int_value();
  const int a = args(3).int_value();
  const int M = args(4).int_value();
  const int N = L/a;
  const int W = coef.columns()/N;

  Matrix f(L,W);  
  
  idgtreal_fac((ltfat_complex*)coef.data(),(ltfat_complex*)gf.data(),L,W,a,M,
	  (double*)f.data());

  return octave_value (f);
}
