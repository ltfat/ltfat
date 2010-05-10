#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_dgtreal_long, args, ,
  "This function calls the C-library\n\
  c=comp_dgtreal_long(f,gf,a,M);\n\
  Yeah.")
{

  const int a = args(2).int_value();
  const int M = args(3).int_value();

  const Matrix f = args(0).matrix_value();
  const Matrix g = args(1).matrix_value();
  
  const int L = f.rows();
  const int W = f.columns();
  const int N = L/a;

  const int M2 = (M/2)+1;
  
  ComplexMatrix cout(M2,N*W);  
  
  dgtreal_long((double*)f.data(),(double*)g.data(),
	    L, W, a, M,(ltfat_complex*)cout.data());
  
  return octave_value (cout);

}
