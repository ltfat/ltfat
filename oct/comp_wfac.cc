#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_wfac, args, ,
  "Computes window factorization.\n\
  Usage: c=comp_wfac(g,a,M);\n\
  Yeah.")
{

  const bool g_is_complex  = args(0).is_complex_type();

  const int a = args(1).int_value();
  const int M = args(2).int_value();

  int h_a, h_m;
   
  if (g_is_complex)
  {
     const ComplexMatrix g = args(0).complex_matrix_value();
     const int L = g.rows();
     
     const int b = L/M;
         
     const int c=gcd(a, M,&h_a, &h_m);
     const int p=a/c;
     const int q=M/c;
     const int d=b/p;
     
     ComplexMatrix gf(p*q,c*d);
     
     wfac((const ltfat_complex*)g.data(),L,a,M,
	  (ltfat_complex*)gf.fortran_vec());
     
     return octave_value (gf);

  }
  else
  {

     const Matrix g = args(0).matrix_value();
     const int L = g.rows();
     const int b = L/M;
         
     const int c=gcd(a, M,&h_a, &h_m);
     const int p=a/c;
     const int q=M/c;
     const int d=b/p;
     
     ComplexMatrix gf(p*q,c*d);
     
     wfac_r(g.data(),L,a,M,(ltfat_complex*)gf.fortran_vec());
     
     return octave_value (gf);



  }
}
