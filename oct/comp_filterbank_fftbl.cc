#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_filterbank_fftbl // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_filterbank_fftbl(F,G,foff,a,realonly)\n Yeah."


#include "ltfat_oct_template_helper.h"


static inline void
fwd_filterbank_fftbl(const Complex *F, const Complex *G[],
                     const ltfatInt L, const ltfatInt Gl[],
                     const ltfatInt W, const double afrac[],
                     const ltfatInt M, const ltfatInt foff[],
                     const int realonly[], Complex *c[])
{
    filterbank_fftbl_d(reinterpret_cast<const fftw_complex *>(F),
                       reinterpret_cast<const fftw_complex **>(G),
                       L, Gl, W, afrac, M, foff, realonly,
                       reinterpret_cast<fftw_complex **>(c));
}

static inline void
fwd_filterbank_fftbl(const FloatComplex *F, const FloatComplex *G[],
                     const ltfatInt L, const ltfatInt Gl[],
                     const ltfatInt W, const double afrac[],
                     const ltfatInt M, const ltfatInt foff[],
                     const int realonly[], FloatComplex *c[])
{

    filterbank_fftbl_s(reinterpret_cast<const fftwf_complex *>(F),
                       reinterpret_cast<const fftwf_complex **>(G),
                       L, Gl, W, afrac, M, foff, realonly,
                       reinterpret_cast<fftwf_complex **>(c));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    // Input data
    MArray<LTFAT_TYPE> F = ltfatOctArray<LTFAT_TYPE>(args(0));
    // Cell aray containing impulse responses
    Cell G = args(1).cell_value();
    Matrix foffDouble = args(2).matrix_value();
    Matrix a = args(3).matrix_value();
    Matrix realonlyDouble = args(4).matrix_value();

    // Input length
    const octave_idx_type L  = F.rows();
    // Number of channels
    const octave_idx_type W  = F.columns();
    // Number of filters
    const octave_idx_type M = G.nelem();


    OCTAVE_LOCAL_BUFFER (double, afrac, M);
    if ( a.columns() > 1)
    {
        for (octave_idx_type m = 0; m < M; m++)
            afrac[m] = a(m) / a(m + M);
    }
    else
    {
        for (octave_idx_type m = 0; m < M; m++)
            afrac[m] = a(m);
    }

    // Allocating temporary arrays
    // Filter lengts
    OCTAVE_LOCAL_BUFFER (ltfatInt, Gl, M);
    // Output subband lengths
    OCTAVE_LOCAL_BUFFER (ltfatInt, foff, M);
    OCTAVE_LOCAL_BUFFER (int, realonly, M);
    // Impulse responses pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, GPtrs, M);
    // Output subbands pointers
    OCTAVE_LOCAL_BUFFER (LTFAT_TYPE*, cPtrs, M);
    //
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, g_elems, M);
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, c_elems, M);

    for (octave_idx_type m = 0; m < M; m++)
    {
        realonly[m] = (int) (realonlyDouble(m) > 1e-3);
        foff[m] = (ltfatInt) foffDouble(m);
        g_elems[m] = ltfatOctArray<LTFAT_TYPE>(G.elem(m));
        GPtrs[m] = g_elems[m].data();
        Gl[m] = (ltfatInt) g_elems[m].numel();
        octave_idx_type outLen = (octave_idx_type) floor( L / afrac[m] + 0.5);
        c_elems[m] = MArray<LTFAT_TYPE>(dim_vector(outLen, W));
        cPtrs[m] = c_elems[m].fortran_vec();
    }

    fwd_filterbank_fftbl(F.data(), GPtrs, L, Gl, W, afrac, M, foff,
                         realonly, cPtrs);

    Cell c(dim_vector(M, 1));
    for (octave_idx_type m = 0; m < M; ++m)
    {
        c.elem(m) = c_elems[m];
    }

    return octave_value(c);
}
