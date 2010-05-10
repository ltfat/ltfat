#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_col2diag, args, ,
  "Computes spreading permutation.\n\
  Usage: cout=comp_col2diag(cin);\n\
  Yeah.")
{

  const bool cin_is_complex  = args(0).is_complex_type();

  if (cin_is_complex)
  {
     const ComplexMatrix cin = args(0).complex_matrix_value();
     const int L = cin.rows();
     
     ComplexMatrix cout(L,L);
     
     col2diag((ltfat_complex*)cin.data(),L,(ltfat_complex*)cout.data());
     
     return octave_value (cout);

  }
  else
  {

     const Matrix cin = args(0).matrix_value();
     const int L = cin.rows();
     
     Matrix cout(L,L);
     
     col2diag_r((double*)cin.data(),L,(double*)cout.data());
     
     return octave_value (cout);



  }
}
