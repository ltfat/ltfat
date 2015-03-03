#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_iatrousfilterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library \n\
                     f=comp_iatrousfilterbank_td(c,g,a,offset) \n\
                     Yeah."

#include "ltfat_oct_template_helper.h"


static inline void
fwd_iatrousfilterbank_td(const Complex *c, const Complex *g[],
                         const ltfatInt L, const ltfatInt gl[],
                         const ltfatInt W, const ltfatInt a[],
                         const ltfatInt offset[], const ltfatInt M,
                         Complex *f, ltfatExtType ext)
{
    iatrousfilterbank_td_cd(reinterpret_cast<const fftw_complex *>(c),
                            reinterpret_cast<const fftw_complex **>(g),
                            L, gl, W, a, offset, M,
                            reinterpret_cast<fftw_complex *>(f),
                            ext);
}

static inline void
fwd_iatrousfilterbank_td(const FloatComplex *c, const FloatComplex *g[],
                         const ltfatInt L, const ltfatInt gl[],
                         const ltfatInt W, const ltfatInt a[],
                         const ltfatInt offset[], const ltfatInt M,
                         FloatComplex *f, ltfatExtType ext)
{
    iatrousfilterbank_td_cs(reinterpret_cast<const fftwf_complex *>(c),
                            reinterpret_cast<const fftwf_complex **>(g),
                            L, gl, W, a, offset, M,
                            reinterpret_cast<fftwf_complex *>(f),
                            ext);
}

static inline void
fwd_iatrousfilterbank_td(const double *c, const double *g[],
                         const ltfatInt L, const ltfatInt gl[],
                         const ltfatInt W, const ltfatInt a[],
                         const ltfatInt offset[], const ltfatInt M,
                         double *f, ltfatExtType ext)
{
    iatrousfilterbank_td_d(c, g, L, gl, W, a, offset, M, f, ext);
}

static inline void
fwd_iatrousfilterbank_td(const float *c, const float *g[],
                         const ltfatInt L, const ltfatInt gl[],
                         const ltfatInt W, const ltfatInt a[],
                         const ltfatInt offset[], const ltfatInt M,
                         float *f, ltfatExtType ext)
{
    iatrousfilterbank_td_s(c, g, L, gl, W, a, offset, M, f, ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    // Input data
    MArray<LTFAT_TYPE> c = ltfatOctArray<LTFAT_TYPE>(args(0));
    // Cell aray containing impulse responses
    MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
    // Subsampling factors
    Matrix aDouble = args(2).matrix_value();
    Matrix offsetDouble = args(3).matrix_value();

    const octave_idx_type L = c.dim1();
    const octave_idx_type M = c.dim2();
    octave_idx_type W = 1;
    if (c.ndims() > 2)
    {
        W = c.dim3();
    }

    const octave_idx_type filtLen = g.rows();

    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, filtLens, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, a, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, offset, M);


    for (octave_idx_type m = 0; m < M; m++)
    {
        filtLens[m] = (ltfatInt) filtLen;
        a[m] = (ltfatInt) aDouble(0);
        offset[m] = (ltfatInt) offsetDouble(m);
        gPtrs[m] = g.data() +  m * filtLen;
    }

    MArray<LTFAT_TYPE> f(dim_vector(L, W));

    fwd_iatrousfilterbank_td(c.data(), gPtrs, L, filtLens, W, a, offset, M,
                             f.fortran_vec(), PER);

    return octave_value(f);
}
