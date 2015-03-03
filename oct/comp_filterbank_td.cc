#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_filterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_filterbank_td(f,g,a,offset,ext)\n Yeah."

#include "ltfat_oct_template_helper.h"

static inline void
fwd_filterbank_td(const Complex *f, const Complex *g[],
                  const ltfatInt L, const ltfatInt gl[],
                  const ltfatInt W, const ltfatInt a[],
                  const ltfatInt offset[], const ltfatInt M,
                  Complex *c[], ltfatExtType ext)
{
    filterbank_td_cd(reinterpret_cast<const fftw_complex *>(f),
                     reinterpret_cast<const fftw_complex **>(g),
                     L, gl, W, a, offset, M,
                     reinterpret_cast<fftw_complex **>(c),
                     ext);
}

static inline void
fwd_filterbank_td(const FloatComplex *f, const FloatComplex *g[],
                  const ltfatInt L, const ltfatInt gl[],
                  const ltfatInt W, const ltfatInt a[],
                  const ltfatInt offset[], const ltfatInt M,
                  FloatComplex *c[], ltfatExtType ext)
{
    filterbank_td_cs(reinterpret_cast<const fftwf_complex *>(f),
                     reinterpret_cast<const fftwf_complex **>(g),
                     L, gl, W, a, offset, M,
                     reinterpret_cast<fftwf_complex **>(c),
                     ext);
}

static inline void
fwd_filterbank_td(const double *f, const double *g[],
                  const ltfatInt L, const ltfatInt gl[],
                  const ltfatInt W, const ltfatInt a[],
                  const ltfatInt offset[], const ltfatInt M,
                  double *c[], ltfatExtType ext)
{
    filterbank_td_d(f, g, L, gl, W, a, offset, M, c, ext);
}

static inline void
fwd_filterbank_td(const float *f, const float *g[],
                  const ltfatInt L, const ltfatInt gl[],
                  const ltfatInt W, const ltfatInt a[],
                  const ltfatInt offset[], const ltfatInt M,
                  float *c[], ltfatExtType ext)
{
    filterbank_td_s(f, g, L, gl, W, a, offset, M, c, ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    // Input data
    MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
    // Cell aray containing impulse responses
    Cell g = args(1).cell_value();
    // Subsampling factors
    Matrix aDouble = args(2).matrix_value();
    // Skips
    Matrix offsetDouble = args(3).matrix_value();
    charMatrix extMat = args(4).char_matrix_value();
    ltfatExtType ext = ltfatExtStringToEnum(extMat.row_as_string(0).c_str());
    // Input length
    const octave_idx_type L  = f.rows();
    // Number of channels
    const octave_idx_type W  = f.columns();
    // Number of filters
    const octave_idx_type M = g.nelem();

    // Allocating temporary arrays
    // Filter lengts
    OCTAVE_LOCAL_BUFFER (ltfatInt, filtLen, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, a, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, offset , M);
    // Impulse responses pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);
    // Output subbands pointers
    OCTAVE_LOCAL_BUFFER (LTFAT_TYPE*, cPtrs, M);
    // Output cell elements array,
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, c_elems, M);
    //
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, g_elems, M);

    for (octave_idx_type m = 0; m < M; m++)
    {
        a[m] = (ltfatInt) aDouble(m);
        offset[m] = (ltfatInt) offsetDouble(m);
        g_elems[m] = ltfatOctArray<LTFAT_TYPE>(g.elem(m));
        gPtrs[m] = g_elems[m].data();
        filtLen[m] = (ltfatInt) g_elems[m].numel();
        octave_idx_type outLen = (octave_idx_type)
                                 filterbank_td_size(L, a[m], filtLen[m],
                                         offset[m], ext);
        c_elems[m] = MArray<LTFAT_TYPE>(dim_vector(outLen, W));
        cPtrs[m] = c_elems[m].fortran_vec();
    }

    fwd_filterbank_td(f.data(), gPtrs, L, filtLen, W, a, offset, M, cPtrs, ext);

    Cell c(dim_vector(M, 1));
    for (octave_idx_type m = 0; m < M; ++m)
    {
        c.elem(m) = c_elems[m];
    }
    return octave_value(c);
}
