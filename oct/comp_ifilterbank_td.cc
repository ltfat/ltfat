#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_ifilterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_ifilterbank_td(c,g,a,Ls,offset,ext);\n Yeah."

#include "ltfat_oct_template_helper.h"


static inline void
fwd_ifilterbank_td(const Complex *c[],  const Complex *g[],
                   const ltfatInt L, const ltfatInt gl[],
                   const ltfatInt W, const ltfatInt a[],
                   const ltfatInt offset[], const ltfatInt M,
                   Complex *f, ltfatExtType ext)
{
    ifilterbank_td_cd(reinterpret_cast<const fftw_complex **>(c),
                      reinterpret_cast<const fftw_complex **>(g),
                      L, gl, W, a, offset, M,
                      reinterpret_cast<fftw_complex *>(f),
                      ext);
}

static inline void
fwd_ifilterbank_td(const FloatComplex *c[],  const FloatComplex *g[],
                   const ltfatInt L, const ltfatInt gl[],
                   const ltfatInt W, const ltfatInt a[],
                   const ltfatInt offset[], const ltfatInt M,
                   FloatComplex *f, ltfatExtType ext)
{
    ifilterbank_td_cs(reinterpret_cast<const fftwf_complex **>(c),
                      reinterpret_cast<const fftwf_complex **>(g),
                      L, gl, W, a, offset, M,
                      reinterpret_cast<fftwf_complex *>(f),
                      ext);
}

static inline void
fwd_ifilterbank_td(const double *c[],  const double *g[],
                   const ltfatInt L, const ltfatInt gl[],
                   const ltfatInt W, const ltfatInt a[],
                   const ltfatInt offset[], const ltfatInt M,
                   double *f, ltfatExtType ext)
{
    ifilterbank_td_d(c, g, L, gl, W, a, offset, M, f, ext);
}

static inline void
fwd_ifilterbank_td(const float *c[],  const float *g[],
                   const ltfatInt L, const ltfatInt gl[],
                   const ltfatInt W, const ltfatInt a[],
                   const ltfatInt offset[], const ltfatInt M,
                   float *f, ltfatExtType ext)
{
    ifilterbank_td_s(c, g, L, gl, W, a, offset, M, f, ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    //DEBUGINFO;
    // Input data
    Cell c = args(0).cell_value();
    // Cell aray containing impulse responses
    Cell g = args(1).cell_value();
    // Subsampling factors
    Matrix aDouble = args(2).matrix_value();
    // Skips
    const octave_idx_type L = args(3).int_value();
    Matrix offsetDouble = args(4).matrix_value();

    charMatrix extMat = args(5).char_matrix_value();
    ltfatExtType ext = ltfatExtStringToEnum(extMat.row_as_string(0).c_str());
    // Number of filters
    const octave_idx_type M = g.nelem();

    // Allocating temporary arrays
    // Filter lengts
    OCTAVE_LOCAL_BUFFER (ltfatInt, filtLen, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, a, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, offset, M);
    // Impulse responses pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);
    // Output subbands pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, cPtrs, M);
    // Input cell elements array,
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, c_elems, M);
    //
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, g_elems, M);

    for (octave_idx_type m = 0; m < M; m++)
    {
        a[m] = (ltfatInt) aDouble(m);
        offset[m] = (ltfatInt) offsetDouble(m);
        g_elems[m] = ltfatOctArray<LTFAT_TYPE>(g.elem(m));
        c_elems[m] = ltfatOctArray<LTFAT_TYPE>(c.elem(m));
        gPtrs[m] = g_elems[m].data();
        cPtrs[m] = c_elems[m].data();
        filtLen[m] = (ltfatInt) g_elems[m].nelem();
    }

    const octave_idx_type W  = c_elems[0].columns();

    MArray<LTFAT_TYPE> f(dim_vector(L, W));

    fwd_ifilterbank_td(cPtrs, gPtrs, L, filtLen, W, a, offset, M,
                       f.fortran_vec(), ext);

    return octave_value(f);
}
