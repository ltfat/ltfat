#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"
#include "math.h"

#define PI 3.1415926535897932384626433832795

DEFUN_DLD (comp_pchirp, args, ,
  "This function calls the C-library\n\
  c=pchirp(L,n);\n")
{
   
   const int L = args(0).int_value();
   const int n = args(1).int_value();
   
   ComplexMatrix g(L,1);
   
   double *gp = (double*)g.fortran_vec();

   const double LL=2.0*L;
   const double Lpone=L+1;
   
   for (int m=0;m<L;m++)
   {
      const double work = PI*fmod(Lpone*n*m*m,LL)/L;
      gp[2*m]   = cos(work);
      gp[2*m+1] = sin(work);
   }


  return octave_value (g);

}
