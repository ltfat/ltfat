#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"


DEFUN_DLD (comp_dgt_long, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_long(f,g,a,M);\n\
  Yeah.")
{

  //const bool f_is_complex  = args(0).is_complex_type();

  const int a = args(2).int_value();
  const int M = args(3).int_value();

  const ComplexMatrix f = args(0).complex_matrix_value();
  const ComplexMatrix g = args(1).complex_matrix_value();
  
  const int L = f.rows();
  const int W = f.columns();
  const int N = L/a;

  dim_vector dims_out(M,N,W);  
  dims_out.chop_trailing_singletons();

  ComplexNDArray cout(dims_out);  
  
  dgt_long((ltfat_complex*)f.data(),(ltfat_complex*)g.data(),
	   L, W, a, M,(ltfat_complex*)cout.data());
  
  return octave_value (cout);


  /*
  if (f_is_complex)
  {


  }
  else
  {

    const Matrix f = args(0).matrix_value();
    ltfat_complex *gf;

    const int L = f.rows();
    const int W = f.columns();
    const int N = L/a;

    gf = (ltfat_complex*)ltfat_malloc(L*sizeof(ltfat_complex));
    
    if (args(1).is_complex_type())
    {
       const ComplexMatrix g = args(1).complex_matrix_value();
       wfac((ltfat_complex*)g.data(),L,a,M,gf);
    }
    else
    {
       const Matrix g = args(1).matrix_value();
       wfac_r((double*)g.data(),L,a,M,gf);
    }
    
    ComplexMatrix cout(M,N*W);  
    
    dgt_fac_r((double*)f.data(),gf,L, W, a, M,
	      (ltfat_complex*)cout.data());
    
    ltfat_free(gf);

    return octave_value (cout);

  }

  */

}
