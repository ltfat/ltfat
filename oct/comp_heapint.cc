#define TYPEDEPARGS 0, 1, 2
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_heapint // change to filename
#define OCTFILEHELP "Computes heapint.\n\
                    Usage: c = comp_heapint(s, itime, ifreq, a, tol);\n\
                    Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer


static inline void fwd_heapint(const double *s, const double *tgrad, const double *fgrad,
                                  const octave_idx_type a, const octave_idx_type M,
								  const octave_idx_type L, const octave_idx_type W,
								  double tol, double *phase)
{
   heapint_d(s,tgrad,fgrad,a,M,L,W,tol,phase);
}

static inline void fwd_heapint(const float *s, const float *tgrad, const float *fgrad,
                               const octave_idx_type a, const octave_idx_type M,
							   const octave_idx_type L, const octave_idx_type W,
							   float tol, float *phase)
{
   heapint_s(s,tgrad,fgrad,a,M,L,W,tol,phase);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 
	 MArray<LTFAT_TYPE> s = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> tgrad = ltfatOctArray<LTFAT_TYPE>(args(1));
	 MArray<LTFAT_TYPE> fgrad = ltfatOctArray<LTFAT_TYPE>(args(2));
	 const octave_idx_type a  = args(3).int_value();
     const double tol   = args(4).double_value();

     const octave_idx_type M = args(0).rows();
     const octave_idx_type N = args(0).columns();
     const octave_idx_type L = N*a;

     MArray<LTFAT_TYPE> phase(dim_vector(M,N)); 
     phase.fill(0);
	 
	 fwd_heapint(s.data(),tgrad.data(),fgrad.data(),a,M,L,1,tol,phase.fortran_vec());
	 
     return octave_value(phase);
}

/*
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
*/