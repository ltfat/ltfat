#include <stdio.h>
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_nonsepdgt_shearola, args, ,
  "This function calls the C-library\n\
  c=comp_nonsepdgt_shearola(f,g,a,M,s0,s1,br,bl);\n")
{

   const ComplexMatrix f = args(0).complex_matrix_value();
   const ComplexMatrix g = args(1).complex_matrix_value();
   const int    a        = args(2).int_value();
   const int    M        = args(3).int_value();
   const int    s0       = args(4).int_value();
   const int    s1       = args(5).int_value();
   const int    br       = args(6).int_value();
   const int    bl       = args(7).int_value();
   

   const int L = f.rows();
   const int W = f.cols();
   const int N = L/a;
   const int gl = g.rows();
  
   dim_vector dims_out(M,N,W);  
   dims_out.chop_trailing_singletons();

   ComplexNDArray cout(dims_out);

   dgt_shearola((const ltfat_complex*)f.fortran_vec(),
		(const ltfat_complex*)g.fortran_vec(),
		L,gl,W,a,M,s0,s1,br,bl,
		(ltfat_complex*)cout.data());
        
   return octave_value (cout);
}
