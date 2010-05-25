#include <octave/oct.h>
#include "config.h"
#include "../dgt_fb_2.c"
#include "../ltfat_time.c"

DEFUN_DLD (mex_dgt_fb_2, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_fb(f,g,a,M,boundary);\n\
  Yeah.")
{

  ComplexMatrix f = args(0).complex_matrix_value();
  ComplexMatrix g = args(1).complex_matrix_value();
  const int a = args(2).int_value();
  const int M = args(3).int_value();
  const int dotime = args(4).int_value();

  const int L  = f.rows();
  const int W  = f.columns();
  const int gl = g.rows();
  const int R  = g.columns();

  const int N = L/a;

  ComplexMatrix cout(M,N*W*R);  
  
  dgt_fb((ltfat_complex*)f.data(),(ltfat_complex*)g.data(),L,gl,W,R,a,M,
	 (ltfat_complex*)cout.data(),dotime);  
  
  return octave_value (cout);
}
