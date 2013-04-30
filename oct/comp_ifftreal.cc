#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_ifftreal // change to filename
#define OCTFILEHELP "This function calls the FFTW3 real FFT\n\
                    f=comp_ifftreal(c,L2);\n\
                    Yeah."


#include "ltfat_oct_template_helper.h"
#include "config.h"
#include "fftw3.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_ifftreal(const Complex *c,
                               const octave_idx_type L,
							   const octave_idx_type W,
							   const octave_idx_type L2,
                               double *f)
{
  fftw_plan p;

  p = fftw_plan_many_dft_c2r(1, &L, W,
			     (double (*)[2])c, NULL,
			     1, L2,
			     (double*)f, NULL,
			     1, L,
			     FFTW_ESTIMATE);
  
  // Real FFT. 
  fftw_execute(p);   
  
  fftw_destroy_plan(p);
}

static inline void fwd_ifftreal(const FloatComplex *c,
                               const octave_idx_type L,
							   const octave_idx_type W,
							   const octave_idx_type L2,
                               float *f)
{
  fftwf_plan p;

  p = fftwf_plan_many_dft_c2r(1, &L, W,
			     (float (*)[2])c, NULL,
			     1, L2,
			     (float*)f, NULL,
			     1, L,
			     FFTW_ESTIMATE);

  // Real FFT. 
  fftwf_execute(p);   
  
  fftwf_destroy_plan(p);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 MArray<LTFAT_TYPE> c = ltfatOctArray<LTFAT_TYPE>(args(0));
  
     const octave_idx_type L2 = c.rows();
     const octave_idx_type W  = c.columns();
     const octave_idx_type L = args(1).int_value();
  
     MArray<LTFAT_REAL> f(dim_vector(L,W)); 
     f.fill(0);
	 
	 fwd_ifftreal(c.data(),L,W,L2,f.fortran_vec());
	 
	 // Scale, because Octave's normalization is different. 
     LTFAT_REAL s  = (LTFAT_REAL) (1.0)/L;
     LTFAT_REAL *fp = f.fortran_vec();
     for (octave_idx_type ii=0; ii<L*W; ii++)
     {
        fp[ii] *=s;
     }
	 
     return octave_value(f);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "fftw3.h"

DEFUN_DLD (comp_ifftreal, args, ,
  "This function calls the FFTW3 real FFT\n\
  f=comp_fftreal(c,L2);\n\
  Yeah.")
{

  fftw_plan p;

  const ComplexMatrix c = args(0).complex_matrix_value();

  const int L2 = c.rows();
  const int W  = c.columns();
  
  const int L = args(1).int_value();
  
  Matrix f(L,W);  
  
  // Prototype:
  //   fftw_plan fftw_plan_many_dft_c2r(int rank, const int *n, int howmany,
  //                                    fftw_complex *in, const int *inembed,
  //                                    int istride, int idist,
  //                                    double *out, const int *onembed,
  //                                    int ostride, int odist,
  //                                    unsigned flags);
  

  // Create plan. Copy data from f to cout. 
  p = fftw_plan_many_dft_c2r(1, &L, W,
			     (fftw_complex*)c.data(), NULL,
			     1, L2,
			     f.fortran_vec(), NULL,
			     1, L,
			     FFTW_ESTIMATE);
  
  // Real IFFT. 
  fftw_execute(p);   

  // Scale, because Octave's normalization is different. 
  double s  = 1.0/L;
  double *fp = f.fortran_vec();
  for (int ii=0; ii<L*W; ii++)
  {
      fp[ii] *=s;
  }
  fftw_destroy_plan(p);   
    
  return octave_value (f);

}
*/