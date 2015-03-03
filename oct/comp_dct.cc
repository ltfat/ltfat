#define TYPEDEPARGS 0
#define REALARGS
#define COMPLEXARGS
#define SINGLEARGS
#define OCTFILENAME comp_dct // change to filename
#define OCTFILEHELP "This function calls FFTW library.\n Yeah."

#include "ltfat_oct_template_helper.h"
#include "config.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void
fwd_dct(const double *f, const octave_idx_type L,
        const octave_idx_type W, double *c,
        const dct_kind kind)
{
    dct_d(f, L, W, c, kind);
}

static inline void
fwd_dct(const float *f, const octave_idx_type L,
        const octave_idx_type W, float *c,
        const dct_kind kind)
{
    dct_s(f, L, W, c, kind);
}

static inline void
fwd_dct(const Complex *f, const octave_idx_type L,
        const octave_idx_type W, Complex *c,
        const dct_kind kind)
{
    dct_cd(reinterpret_cast<const fftw_complex *>(f),
           L, W,
           reinterpret_cast<fftw_complex *>(c),
           kind);
}

static inline void
fwd_dct(const FloatComplex *f, const octave_idx_type L,
        const octave_idx_type W, FloatComplex* c,
        const dct_kind kind)
{
    dct_cs(reinterpret_cast<const fftwf_complex *>(f),
           L, W,
           reinterpret_cast<fftwf_complex *>(c),
           kind);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    DEBUGINFO;
    dct_kind kind = DCTI;
    const octave_idx_type type = args(1).int_value();
    MArray<LTFAT_TYPE> f = MArray<LTFAT_TYPE>(ltfatOctArray<LTFAT_TYPE>(args(0)));
    const octave_idx_type L  = f.rows();
    const octave_idx_type W  = f.columns();
    MArray<LTFAT_TYPE> c(dim_vector(L, W));

    switch (type)
    {
    case 1:
        kind = DCTI;
        break;
    case 2:
        kind = DCTII;
        break;
    case 3:
        kind = DCTIII;
        break;
    case 4:
        kind = DCTIV;
        break;
    default:
        error("Unknown type.");
    }

    fwd_dct(f.data(), L, W, c.fortran_vec(), kind);

    return octave_value(c);
}


