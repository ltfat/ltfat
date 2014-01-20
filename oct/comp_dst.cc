#define TYPEDEPARGS 0
#define REALARGS
#define COMPLEXARGS
#define SINGLEARGS
#define OCTFILENAME comp_dst // change to filename
#define OCTFILEHELP "This function calls FFTW library.\n Yeah."

#include "ltfat_oct_template_helper.h"
#include "config.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_dst(const double *f,
                           const octave_idx_type L,
                           const octave_idx_type W,
                           const dst_kind kind,
                           double *c)
{
    dst_d(f,L,W,kind,c);
}

static inline void fwd_dst(const float *f,
                           const octave_idx_type L,
                           const octave_idx_type W,
                           const dst_kind kind,
                           float *c)
{
    dst_s(f,L,W,kind,c);
}

static inline void fwd_dst(const Complex *f,
                           const octave_idx_type L,
                           const octave_idx_type W,
                           const dst_kind kind,
                           Complex *c)
{
    dst_cd(reinterpret_cast<const double _Complex *>(f),
           L,W,kind,
           reinterpret_cast<double _Complex *>(c));
}

static inline void fwd_dst(const FloatComplex *f,
                           const octave_idx_type L,
                           const octave_idx_type W,
                           const dst_kind kind,
                           FloatComplex* c)
{
    dst_cs(reinterpret_cast<const float _Complex *>(f),
           L,W,kind,
           reinterpret_cast<float _Complex *>(c));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    DEBUGINFO;
    dst_kind kind = DSTI;
    const octave_idx_type type = args(1).int_value();
    MArray<LTFAT_TYPE> f = MArray<LTFAT_TYPE>(ltfatOctArray<LTFAT_TYPE>(args(0)));
    const octave_idx_type L  = f.rows();
    const octave_idx_type W  = f.columns();
    MArray<LTFAT_TYPE> c(dim_vector(L,W));

    switch(type)
    {
    case 1:
        kind = DSTI;
        break;
    case 2:
        kind = DSTII;
        break;
    case 3:
        kind = DSTIII;
        break;
    case 4:
        kind = DSTIV;
        break;
    default:
        error("Unknown type.");
    }

    fwd_dst(f.data(),L,W,kind,c.fortran_vec());

    return octave_value(c);
}



