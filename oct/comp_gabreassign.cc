#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_gabreassign, args, ,
  "Computes spreading permutation.\n\
  Usage: sr=comp_gabreassign(s,itime,ifreq,a);\n\
  Yeah.")
{

     const Matrix s     = args(0).matrix_value();
     const Matrix tgrad = args(1).matrix_value();
     const Matrix fgrad = args(2).matrix_value();
     const int a  = args(3).int_value();
     const int M  = s.rows();
     const int N  = s.columns();
     const int L  = N*a;

     Matrix sr(M,N);
     
     gabreassign((double*)s.data(),(double*)tgrad.data(),(double*)fgrad.data(),
		 L,1,a,M,(double*)sr.data());
     
     return octave_value (sr);

}
