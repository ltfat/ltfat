#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_dgt_fb, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_fb(f,g,a,M,boundary);\n\
  Yeah.")
{

  const bool use_complex  = args(0).is_complex_type() || args(1).is_complex_type();

  const int a = args(2).int_value();
  const int M = args(3).int_value();


  if (use_complex)
  {  
     const ComplexMatrix f = args(0).complex_matrix_value();
     const ComplexMatrix g = args(1).complex_matrix_value();
     
     const int L  = f.rows();
     const int W  = f.columns();
     const int gl = g.rows();
     const int N = L/a;

     dim_vector dims_out(M,N,W);  
     dims_out.chop_trailing_singletons();
     
     ComplexNDArray cout(dims_out);  
     
     dgt_fb((ltfat_complex*)f.data(),(ltfat_complex*)g.data(),L,gl,W,a,M,
     (ltfat_complex*)cout.data());  

     return octave_value (cout);
  
  }
  else
  {

     const Matrix f = args(0).matrix_value();
     const Matrix g = args(1).matrix_value();
     
     const int L  = f.rows();
     const int W  = f.columns();
     const int gl = g.rows();
     const int N = L/a;
     
     dim_vector dims_out(M,N,W);  
     dims_out.chop_trailing_singletons();

     ComplexNDArray cout(dims_out);  
     
     dgt_fb_r((double*)f.data(),(double*)g.data(),L,gl,W,a,M,
     (ltfat_complex*)cout.data());  

     return octave_value (cout);


  }


}
