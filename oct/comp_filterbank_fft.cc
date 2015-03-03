#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_filterbank_fft // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_filterbank_fft(F,G,a)\n Yeah."

#include "ltfat_oct_template_helper.h"

static inline void
fwd_filterbank_fft(const Complex *F, const Complex *G[],
                   const ltfatInt L, const ltfatInt W,
                   const ltfatInt a[], const ltfatInt M,
                   Complex *c[])
{
    filterbank_fft_d(reinterpret_cast<const fftw_complex *>(F),
                     reinterpret_cast<const fftw_complex **>(G),
                     L, W, a, M,
                     reinterpret_cast<fftw_complex **>(c));
}

static inline void
fwd_filterbank_fft(const FloatComplex *F, const FloatComplex *G[],
                   const ltfatInt L, const ltfatInt W,
                   const ltfatInt a[], const ltfatInt M,
                   FloatComplex *c[])
{
    filterbank_fft_s(reinterpret_cast<const fftwf_complex *>(F),
                     reinterpret_cast<const fftwf_complex **>(G),
                     L, W, a, M,
                     reinterpret_cast<fftwf_complex **>(c));
}


template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    // Input data
    MArray<LTFAT_TYPE> F = ltfatOctArray<LTFAT_TYPE>(args(0));
    // Cell aray containing impulse responses
    Cell G = args(1).cell_value();
    // Subsampling factors
    Matrix aDouble = args(2).matrix_value();

    // Input length
    const octave_idx_type L  = F.rows();
    // Number of channels
    const octave_idx_type W  = F.columns();
    // Number of filters
    const octave_idx_type M = G.nelem();

    // Allocating temporary arrays
    // Output subband lengths
    OCTAVE_LOCAL_BUFFER (ltfatInt, a, M);
    // Impulse responses pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, GPtrs, M);
    // Output subbands pointers
    OCTAVE_LOCAL_BUFFER (LTFAT_TYPE*, cPtrs, M);
    // Output cell elements array,
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, gElems, M);
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, c_elems, M);

    for (octave_idx_type m = 0; m < M; m++)
    {
        a[m] = (ltfatInt)aDouble(m);
        gElems[m] = ltfatOctArray<LTFAT_TYPE>(G.elem(m));
        GPtrs[m] = gElems[m].data();
        octave_idx_type outLen = (octave_idx_type) ceil( L / aDouble(m) );
        c_elems[m] = MArray<LTFAT_TYPE>(dim_vector(outLen, W));
        cPtrs[m] = c_elems[m].fortran_vec();
    }

    fwd_filterbank_fft(F.data(), GPtrs, L, W, a, M, cPtrs);

    Cell c(dim_vector(M, 1));
    for (octave_idx_type m = 0; m < M; ++m)
    {
        c.elem(m) = c_elems[m];
    }

    return octave_value(c);
}
