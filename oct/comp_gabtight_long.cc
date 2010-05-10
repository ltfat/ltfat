#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_gabtight_long, args, ,
  "This function calls the C-library\n\
  gd=comp_gabtight_long(g,a,M);\n\
  Yeah.")
{

  if (args(0).is_complex_type())
  {

     const ComplexMatrix g = args(0).complex_matrix_value();
     const int L = g.rows();
     const int a = args(1).int_value();
     const int M = args(2).int_value();

     ComplexMatrix gd(L,1);
     
     gabtight_long((const ltfat_complex*)g.fortran_vec(),
		  L,a,M,
		  (ltfat_complex*)gd.data());

     return octave_value (gd);

  }
  else
  {

     const Matrix g = args(0).matrix_value();
     const int L = g.rows();
     const int a = args(1).int_value();
     const int M = args(2).int_value();

     Matrix gd(L,1);
  
     gabtightreal_long((double*)g.fortran_vec(),
		  L,a,M,
		  (double*)gd.fortran_vec());

     return octave_value (gd);

  }
}
