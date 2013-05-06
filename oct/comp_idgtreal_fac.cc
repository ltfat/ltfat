#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_idgtreal_fac // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                    c=comp_idgtreal_fac(c,gf,L,W,a,M);\n\
                    Yeah.\n"


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_idgtreal_fac(const Complex *coef, const Complex *gf,
                               const octave_idx_type L, const octave_idx_type W,
							   const octave_idx_type a, const octave_idx_type M,
							   double *f)
{
   idgtreal_fac_d(reinterpret_cast<const double (*)[2]>(coef),
              reinterpret_cast<const double (*)[2]>(gf),
              L,W,a,M,f);
}

static inline void fwd_idgtreal_fac(const FloatComplex *coef, const FloatComplex *gf,
                               const octave_idx_type L, const octave_idx_type W, 
							   const octave_idx_type a, const octave_idx_type M,
							   float *f)
{
   idgtreal_fac_s(reinterpret_cast<const float (*)[2]>(coef),
              reinterpret_cast<const float (*)[2]>(gf),
              L,W,a,M,f);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 
	 MArray<LTFAT_TYPE> coef = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> gf = ltfatOctArray<LTFAT_TYPE>(args(1));
	 const octave_idx_type L = args(2).int_value();
     const octave_idx_type a = args(3).int_value();
     const octave_idx_type M = args(4).int_value();
     const octave_idx_type N = L/a;
     const octave_idx_type W = coef.columns()/N;

     MArray<LTFAT_REAL> f(dim_vector(L,W)); 
     f.fill(0);
	 
	 fwd_idgtreal_fac(coef.data(),gf.data(),L,W,a,M,f.fortran_vec());
	 
     return octave_value(f);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_idgtreal_fac, args, ,
  "This function calls the C-library\n\
  c=comp_idgt_fac(c,gf,L,a,M);\n\
  Yeah.\n")
{

  const ComplexMatrix coef = args(0).complex_matrix_value();
  const ComplexMatrix gf = args(1).complex_matrix_value();
  const int L = args(2).int_value();
  const int a = args(3).int_value();
  const int M = args(4).int_value();
  const int N = L/a;
  const int W = coef.columns()/N;

  Matrix f(L,W);  
  
  idgtreal_fac((ltfat_complex*)coef.data(),(ltfat_complex*)gf.data(),L,W,a,M,
	  (double*)f.data());

  return octave_value (f);
}
*/