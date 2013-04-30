#define TYPEDEPARGS 0
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_fftreal // change to filename
#define OCTFILEHELP "This function calls the FFTW3 real FFT\n\
                     c=comp_fftreal(f);\n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
#include "config.h"
#include "fftw3.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_fftreal(const double *f,
                               const octave_idx_type L,
							   const octave_idx_type W,
							   const octave_idx_type L2,
                               Complex *cout)
{
  fftw_plan p;
  // Create plan. Copy data from f to cout. 
  p = fftw_plan_many_dft_r2c(1, &L, W,
			     (double*)f, NULL,
			     1, L,
			     (double (*)[2])cout, NULL,
			     1, L2,
			     FFTW_OPTITYPE);

  // Real FFT. 
  fftw_execute(p);   
  
  fftw_destroy_plan(p);
}

static inline void fwd_fftreal(const float *f,
                               const octave_idx_type L,
							   const octave_idx_type W,
							   const octave_idx_type L2,
                               FloatComplex *cout)
{
  fftwf_plan p;
  // Create plan. Copy data from f to cout. 
  p = fftwf_plan_many_dft_r2c(1, &L, W,
			     (float*)f, NULL,
			     1, L,
			     (float (*)[2])cout, NULL,
			     1, L2,
			     FFTW_OPTITYPE);

  // Real FFT. 
  fftwf_execute(p);   
  
  fftwf_destroy_plan(p);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
  
     const octave_idx_type L = f.rows();
     const octave_idx_type W = f.columns();
     const octave_idx_type L2 = (L/2)+1;
  
     MArray<LTFAT_COMPLEX> cout(dim_vector(L2,W)); 
     cout.fill(0);
	 
	 fwd_fftreal(f.data(),L,W,L2,cout.fortran_vec());
	 
     return octave_value(cout);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "fftw3.h"

DEFUN_DLD (comp_fftreal, args, ,
  "This function calls the FFTW3 real FFT\n\
  c=comp_fftreal(f);\n\
  Yeah.")
{

  fftw_plan p;

  const Matrix f = args(0).matrix_value();

  const int L = f.rows();
  const int W = f.columns();
  
  const int L2 = (L/2)+1;
  
  ComplexMatrix cout(L2,W);  
  
  // Create plan. Copy data from f to cout. 
  p = fftw_plan_many_dft_r2c(1, &L, W,
			     (double*)f.data(), NULL,
			     1, L,
			     (fftw_complex*)cout.data(), NULL,
			     1, L2,
			     FFTW_OPTITYPE);
  
  // Real FFT. 
  fftw_execute(p);   
  
  fftw_destroy_plan(p);   
    
  return octave_value (cout);

}
*/