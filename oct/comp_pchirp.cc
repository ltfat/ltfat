#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_pchirp, args, ,
  "This function calls the C-library\n\
  c=pchirp(L,n);\n")
{
   
   const int L = args(0).int_value();
   const int n = args(1).int_value();
   
   ComplexMatrix g(L,1);
   
   pchirp(L, n, (ltfat_complex*)g.fortran_vec());

   return octave_value (g);

}
