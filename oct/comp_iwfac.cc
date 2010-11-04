#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_iwfac, args, ,
  "Computes inverse window factorization.\n\
  Usage: c=comp_iwfac(gf,L,a,M);\n\
  Yeah.")
{

  const ComplexMatrix gf = args(0).complex_matrix_value();
  const int L = args(1).int_value();
  const int a = args(2).int_value();
  const int M = args(3).int_value();

  ComplexMatrix g(L,1);
  
  iwfac((const ltfat_complex*)gf.data(),L,a,M,
	(ltfat_complex*)g.fortran_vec());

  return octave_value (g);
}
