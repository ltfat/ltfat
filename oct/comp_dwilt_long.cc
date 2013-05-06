#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_dwilt_long // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     coef=comp_dwilt_long(f,g,M,L);\n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer
/*
  dgtreal_ola forwarders
*/

static inline void fwd_dwilt_long(const Complex *f, const Complex *g,
                                  const octave_idx_type L, 
							      const octave_idx_type W, 
							      const octave_idx_type M, 
								  Complex *cout)
{
   dwilt_long_d(reinterpret_cast<const double (*)[2]>(f),
                reinterpret_cast<const double (*)[2]>(g),
				L,W,M,
				reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_dwilt_long(const FloatComplex *f, 
                                  const FloatComplex *g,
                                  const octave_idx_type L, 
							      const octave_idx_type W, 
							      const octave_idx_type M, 
								  FloatComplex *cout)
{
   dwilt_long_s(reinterpret_cast<const float (*)[2]>(f),
                reinterpret_cast<const float (*)[2]>(g),
				L,W,M,
				reinterpret_cast<float (*)[2]>(cout));
}

static inline void fwd_dwilt_long(const double *f, const double *g,
                                  const octave_idx_type L, 
							      const octave_idx_type W, 
							      const octave_idx_type M, 
								  double *cout)
{
   dwiltreal_long_d(f,g,L,W,M,cout);
}

static inline void fwd_dwilt_long(const float *f, const float *g,
                                  const octave_idx_type L, 
							      const octave_idx_type W, 
							      const octave_idx_type M, 
								  float *cout)
{
   dwiltreal_long_s(f,g,L,W,M,cout);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
     const octave_idx_type M = args(2).int_value();
     const octave_idx_type L = args(3).int_value();
     const octave_idx_type N = L/M;

	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
  
     const int W = f.columns();

     dim_vector dims_out(2*M,N/2,W);  
     dims_out.chop_trailing_singletons();

     MArray<LTFAT_TYPE> cout(dims_out); 
     cout.fill(0);
	 
	 fwd_dwilt_long(f.data(),g.data(),L, W, M,cout.fortran_vec());
	 
     return octave_value(cout);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"


DEFUN_DLD (comp_dwilt_long, args, ,
  "This function calls the C-library\n\
  coef=comp_dwilt_long(f,g,M,L);\n\
  Yeah.")
{
  const int M = args(2).int_value();
  const int L = args(3).int_value();
  const int N = L/M;

  if (args(0).is_complex_type())
  {

    const ComplexMatrix f = args(0).complex_matrix_value();
    const ComplexMatrix g = args(1).complex_matrix_value();
    
    const int W = f.columns();

    dim_vector dims_out(2*M,N/2,W);  
    dims_out.chop_trailing_singletons();
    
    ComplexNDArray cout(dims_out);  
    
    dwilt_long((ltfat_complex*)f.data(),(ltfat_complex*)g.data(), L, W, M,
	       (ltfat_complex*)cout.data());

    return octave_value (cout);

  }
  else
  {

    const Matrix f = args(0).matrix_value();
    const Matrix g = args(1).matrix_value();
    
    const int W = f.columns();

    dim_vector dims_out(2*M,N/2,W);  
    dims_out.chop_trailing_singletons();
    
    NDArray cout(dims_out);  
    
    dwiltreal_long((double*)f.data(),(double*)g.data(), L, W, M,
		  (double*)cout.data());
	
    return octave_value (cout);

  }

}
*/