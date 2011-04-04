#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_pfilt_fir, args, ,
  "This function calls the C-library\n\
  c=comp_pfilt_fir(f,g,a);\n\
  Yeah.")
{

  const Matrix f = args(0).matrix_value();
  const Matrix g = args(1).matrix_value();
  const int a = args(2).int_value();
     
  const int L  = f.rows();
  const int W  = f.columns();
  const int gl = g.rows();
  const int N = L/a;
  
  Matrix cout(N,W);
  
  pfilt_fir_rr(f.data(),g.data(),L,gl,W,a,cout.fortran_vec());  
  
  return octave_value (cout);

}
