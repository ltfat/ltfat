#include <octave/oct.h>
#include "config.h"
#include "../dgt_fac_1.c"
#include "../gcd.c"
#include "../ltfat_time.c"

DEFUN_DLD (mex_dgt_fac_1, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_fac(f,gf,a,M,dotime);\n\
  Yeah.")
{

  ComplexMatrix f = args(0).complex_matrix_value();
  ComplexMatrix gf = args(1).complex_matrix_value();
  const int a = args(2).int_value();
  const int M = args(3).int_value();
  const int dotime = args(4).int_value();

  const int L = f.rows();
  const int W = f.columns();
  const int R = gf.rows()*gf.columns()/L;

  const int N = L/a;

  ComplexMatrix cout(M,N*W*R);  
  
  dgt_fac((ltfat_complex*)f.data(),(ltfat_complex*)gf.data(),L,W,R,a,M,
	  (ltfat_complex*)cout.data(),dotime);  

  return octave_value (cout);
}
