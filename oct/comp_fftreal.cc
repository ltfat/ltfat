#define TYPEDEPARGS 0
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_fftreal // change to filename
#define OCTFILEHELP "This function calls the FFTW3 real FFT\n\
                     c=comp_fftreal(f);\n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_fftreal(const double *f,
                               const octave_idx_type L,
                               const octave_idx_type W,
                               Complex *cout)
{
    fftreal_d(const_cast<double*>(f),
              L, W,
              reinterpret_cast<fftw_complex*>(cout));
}

static inline void fwd_fftreal(const float *f,
                               const octave_idx_type L,
                               const octave_idx_type W,
                               FloatComplex *cout)
{
    fftreal_s(const_cast<float*>(f),
              L, W,
              reinterpret_cast<fftwf_complex*>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    DEBUGINFO;
    MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));

    const octave_idx_type L = f.rows();
    const octave_idx_type W = f.columns();
    const octave_idx_type L2 = (L / 2) + 1;

    MArray<LTFAT_COMPLEX> cout(dim_vector(L2, W));

    fwd_fftreal(f.data(), L, W, cout.fortran_vec());

    return octave_value(cout);
}

