#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_heapint, args, ,
  "Computes heapint.\n\
  Usage: c = comp_heapint(s, itime, ifreq, a, tol);\n\
  Yeah.")
{

  const Matrix s     = args(0).matrix_value();
  const Matrix tgrad = args(1).matrix_value();
  const Matrix fgrad = args(2).matrix_value();
  const int a        = args(3).int_value();
  const double tol   = args(4).double_value();

  const int M = args(0).rows();
  const int N = args(0).columns();
  const int L = N*a;

  Matrix phase(M,N);
  
  heapint((double*)s.data(),
	  (double*)tgrad.data(),
	  (double*)fgrad.data(),a,M,L,1,tol,
	  (double*)phase.data());
  
  return octave_value (phase);
}
