#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_dgt_fb // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=comp_dgt_fb(f,g,a,M,boundary);\n Yeah."

#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer
/*
  dgt_fb forwarders
*/

static inline void fwd_dgt_fb(const Complex *f, const Complex *g,
                              const octave_idx_type L, const octave_idx_type gl,
		                      const octave_idx_type W, const octave_idx_type a,
							  const octave_idx_type M, Complex *cout)
{
   dgt_fb_cd(reinterpret_cast<const double (*)[2]>(f),
             reinterpret_cast<const double (*)[2]>(g),
             L,gl,W,a,M,
			 reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_dgt_fb(const FloatComplex *f, const FloatComplex *g,
                              const octave_idx_type L, const octave_idx_type gl,
		                      const octave_idx_type W, const octave_idx_type a,
							  const octave_idx_type M, FloatComplex *cout)
{
   dgt_fb_cs(reinterpret_cast<const float (*)[2]>(f),
             reinterpret_cast<const float (*)[2]>(g),
             L,gl,W,a,M,
			 reinterpret_cast<float (*)[2]>(cout));
}

static inline void fwd_dgt_fb(const double *f, const double *g,
                              const octave_idx_type L, const octave_idx_type gl,
		                      const octave_idx_type W, const octave_idx_type a,
							  const octave_idx_type M, Complex *cout)
{
   dgt_fb_d(reinterpret_cast<const double*>(f),
            reinterpret_cast<const double*>(g),
            L,gl,W,a,M,
			reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_dgt_fb(const float *f, const float *g,
                              const octave_idx_type L, const octave_idx_type gl,
		                      const octave_idx_type W, const octave_idx_type a,
							  const octave_idx_type M, FloatComplex *cout)
{
   dgt_fb_s(reinterpret_cast<const float*>(f),
            reinterpret_cast<const float*>(g),
            L,gl,W,a,M,
			reinterpret_cast<float (*)[2]>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 const octave_idx_type a = args(2).int_value();
     const octave_idx_type M = args(3).int_value();
	 
	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
	 const octave_idx_type L  = f.rows();
     const octave_idx_type W  = f.columns();
     const octave_idx_type gl = g.rows();
     const octave_idx_type N = L/a;
	 
	 dim_vector dims_out(M,N,W);  
     dims_out.chop_trailing_singletons();

     MArray<LTFAT_COMPLEX> cout(dims_out); 
     cout.fill(0);
	 
	 fwd_dgt_fb(f.data(),g.data(),L,gl,W,a,M,cout.fortran_vec());
	 
     return octave_value(cout);
}
/* OLD CODE
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_dgt_fb, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_fb(f,g,a,M,boundary);\n\
  Yeah.")
{

  const bool use_complex  = args(0).is_complex_type() || args(1).is_complex_type();

  const int a = args(2).int_value();
  const int M = args(3).int_value();


  if (use_complex)
  {  
     const ComplexMatrix f = args(0).complex_matrix_value();
     const ComplexMatrix g = args(1).complex_matrix_value();
     
     const int L  = f.rows();
     const int W  = f.columns();
     const int gl = g.rows();
     const int N = L/a;

     dim_vector dims_out(M,N,W);  
     dims_out.chop_trailing_singletons();
     
     ComplexNDArray cout(dims_out);  
     
     dgt_fb((ltfat_complex*)f.data(),(ltfat_complex*)g.data(),L,gl,W,a,M,
     (ltfat_complex*)cout.data());  

     return octave_value (cout);
  
  }
  else
  {

     const Matrix f = args(0).matrix_value();
     const Matrix g = args(1).matrix_value();
     
     const int L  = f.rows();
     const int W  = f.columns();
     const int gl = g.rows();
     const int N = L/a;
     
     dim_vector dims_out(M,N,W);  
     dims_out.chop_trailing_singletons();

     ComplexNDArray cout(dims_out);  
     
     dgt_fb_r((double*)f.data(),(double*)g.data(),L,gl,W,a,M,
     (ltfat_complex*)cout.data());  

     return octave_value (cout);


  }


}
*/