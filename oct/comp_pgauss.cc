#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_pgauss, args, ,
  "This function calls the C-library\n\
  c=comp_pgauss(L,wf,center);\n")
{

  const int    L      = args(0).int_value();
  const double w      = args(1).double_value();
  const double center = args(2).double_value();

  Matrix g(L,1);

  pgauss(L, w, center, g.fortran_vec());

  return octave_value (g);
}
