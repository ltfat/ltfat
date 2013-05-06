#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_iwfac // change to filename
#define OCTFILEHELP "Computes inverse window factorization.\n\
                    Usage: c=comp_iwfac(gf,L,a,M);\n\
                    Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_iwfac( const Complex *gf,
                               const octave_idx_type L, const octave_idx_type R,
							   const octave_idx_type a, const octave_idx_type M,
							   Complex *g)
{
   iwfac_d(reinterpret_cast<const double (*)[2]>(gf),
           L,R,a,M,
		   reinterpret_cast<double (*)[2]>(g));
}

static inline void fwd_iwfac( const FloatComplex *gf,
                               const octave_idx_type L, const octave_idx_type R, 
							   const octave_idx_type a, const octave_idx_type M,
							   FloatComplex *g)
{
   iwfac_s(reinterpret_cast<const float (*)[2]>(gf),
           L,R,a,M,
		   reinterpret_cast<float (*)[2]>(g));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 
	 MArray<LTFAT_TYPE> gf = ltfatOctArray<LTFAT_TYPE>(args(0));
	 const octave_idx_type L = args(1).int_value();
     const octave_idx_type R = gf.numel()/L;
     const octave_idx_type a = args(2).int_value();
     const octave_idx_type M = args(3).int_value();

     MArray<LTFAT_COMPLEX> g(dim_vector(L,1)); 
     g.fill(0);
	 
	 fwd_iwfac(gf.data(),L,R,a,M,g.fortran_vec());
	 
     return octave_value(g);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_iwfac, args, ,
  "Computes inverse window factorization.\n\
  Usage: c=comp_iwfac(gf,L,a,M);\n\
  Yeah.")
{

  const ComplexMatrix gf = args(0).complex_matrix_value();
  const int L = args(1).int_value();
  const int R = gf.numel()/L;
  const int a = args(2).int_value();
  const int M = args(3).int_value();

  ComplexMatrix g(L,1);
  
  iwfac((const ltfat_complex*)gf.data(),L,R,a,M,
	(ltfat_complex*)g.fortran_vec());

  return octave_value (g);
}
*/