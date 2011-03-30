#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_idgt_fb, args, ,
  "This function calls the C-library\n\
  c=comp_idgt_fb(coef,g,L,a,M);\n\
  Yeah.")
{   
   const ComplexMatrix coef = args(0).complex_matrix_value();

   const int L = args(2).int_value();
   const int a = args(3).int_value();
   const int M = args(4).int_value();
   const int N = L/a;

   const int W = coef.rows()*coef.columns()/(M*N);
   
   ComplexMatrix f(L,W);  
   
   if (args(1).is_complex_type())
   {
      const ComplexMatrix g = args(1).complex_matrix_value();
      const int gl = g.rows();
      
      idgt_fb((const ltfat_complex*)coef.data(),(const ltfat_complex*)g.data(),
	      L,gl,W,a,M,
	      (ltfat_complex*)f.fortran_vec());        
   }
   else
   {
      const Matrix g = args(1).matrix_value();
      const int gl = g.rows();
      
      idgt_fb_r((const ltfat_complex*)coef.data(),g.data(),
		L,gl,W,a,M,
		(ltfat_complex*)f.fortran_vec());      
   }

   return octave_value (f);               
}
