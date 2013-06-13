#define TYPEDEPARGS 0
#define SINGLEARGS
#define MATCHEDARGS 1
#define COMPLEXARGS
#define OCTFILENAME comp_idgtreal_fb // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                    c=comp_idgtereal_fb(c,g,L,a,M);\n\
                    Yeah.\n"


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_idgtreal_fb(const Complex *coef, const double *gf,
                               const octave_idx_type L,
                               const octave_idx_type gl, const octave_idx_type W,
							   const octave_idx_type a, const octave_idx_type M,
							   double *f)
{
   idgtreal_fb_d(reinterpret_cast<const double (*)[2]>(coef),
                 gf,L,gl,W,a,M,f);
}

static inline void fwd_idgtreal_fb(const FloatComplex *coef, const float *gf,
                               const octave_idx_type L,
                               const octave_idx_type gl, const octave_idx_type W, 
							   const octave_idx_type a, const octave_idx_type M,
							   float *f)
{
   idgtreal_fb_s(reinterpret_cast<const float (*)[2]>(coef),
                 gf,L,gl,W,a,M,f);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 
	 MArray<LTFAT_TYPE> coef = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_REAL> gf = ltfatOctArray<LTFAT_REAL>(args(1));
	 const octave_idx_type L = args(2).int_value();
     const octave_idx_type a = args(3).int_value();
     const octave_idx_type M = args(4).int_value();
     const octave_idx_type N = L/a;
     const octave_idx_type W = coef.columns()/N;
     const octave_idx_type gl = gf.rows();

     MArray<LTFAT_REAL> f(dim_vector(L,W)); 
     f.fill(0);
	 
	 fwd_idgtreal_fb(coef.data(),gf.data(),L,gl,W,a,M,f.fortran_vec());
	 
     return octave_value(f);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_idgtreal_fb, args, ,
  "This function calls the C-library\n\
  c=comp_idgt_fb(c,g,L,a,M);\n\
  Yeah.\n")
{

  const ComplexMatrix coef = args(0).complex_matrix_value();
  const Matrix g = args(1).matrix_value();
  const int L = args(2).int_value();
  const int a = args(3).int_value();
  const int M = args(4).int_value();
  const int N = L/a;
  const int W = coef.columns()/N;
  const int gl = g.rows();

  Matrix f(L,W);  
  
  idgtreal_fb((ltfat_complex*)coef.data(),g.data(),L,gl,W,a,M,
	  f.fortran_vec());

  return octave_value (f);
}
*/