#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_dgt_ola // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_dgt_ola(f,g,a,M,bl);\n Yeah."


#include "ltfat_oct_template_helper.h"
/*
  dgt_ola forwarders
*/

static inline void
fwd_dgt_ola(const Complex *f, const Complex *g, const octave_idx_type L,
            const octave_idx_type gl, const octave_idx_type W,
            const octave_idx_type a, const octave_idx_type M,
            const octave_idx_type bl, const octave_idx_type ptype, Complex *cout)
{
    dgt_ola_d(reinterpret_cast<const fftw_complex*>(f),
              reinterpret_cast<const fftw_complex*>(g),
              L, gl, W, a, M, bl,
              static_cast<const dgt_phasetype>(ptype),
              reinterpret_cast<fftw_complex*>(cout));
}

static inline void
fwd_dgt_ola(const FloatComplex *f, const FloatComplex *g,
            const octave_idx_type L, const octave_idx_type gl,
            const octave_idx_type W, const octave_idx_type a,
            const octave_idx_type M, const octave_idx_type bl,
            const octave_idx_type ptype, FloatComplex *cout)
{
    dgt_ola_s(reinterpret_cast<const fftwf_complex*>(f),
              reinterpret_cast<const fftwf_complex*>(g),
              L, gl, W, a, M, bl,
              static_cast<const dgt_phasetype>(ptype),
              reinterpret_cast<fftwf_complex*>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list
octFunction(const octave_value_list& args, int nargout)
{
    DEBUGINFO;
    const octave_idx_type a  = args(2).int_value();
    const octave_idx_type M  = args(3).int_value();
    const octave_idx_type bl = args(4).int_value();
    const octave_idx_type ptype = args(4).int_value();

    MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
    MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));

    const octave_idx_type L  = f.rows();
    const octave_idx_type W  = f.columns();
    const octave_idx_type gl = g.rows();
    const octave_idx_type N = L / a;

    dim_vector dims_out(M, N, W);
    dims_out.chop_trailing_singletons();

    MArray<LTFAT_TYPE> cout(dims_out);

    fwd_dgt_ola(f.data(), g.data(), L, gl, W, a, M, bl, ptype,
                cout.fortran_vec());

    return octave_value(cout);
}
