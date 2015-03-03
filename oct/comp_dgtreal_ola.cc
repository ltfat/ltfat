#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_dgtreal_ola // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_dgtreal_ola(f,g,a,M,bl);\n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer
/*
  dgtreal_ola forwarders
*/

static inline void
fwd_dgtreal_ola(const double *f, const double *g,
                const octave_idx_type L, const octave_idx_type gl,
                const octave_idx_type W, const octave_idx_type a,
                const octave_idx_type M, const octave_idx_type bl,
                const octave_idx_type ptype, Complex *cout)
{
    dgtreal_ola_d(f, g, L, gl, W, a, M, bl,
                  static_cast<dgt_phasetype>(ptype),
                  reinterpret_cast<fftw_complex*>(cout));
}

static inline void
fwd_dgtreal_ola(const float *f, const float *g,
                const octave_idx_type L, const octave_idx_type gl,
                const octave_idx_type W, const octave_idx_type a,
                const octave_idx_type M, const octave_idx_type bl,
                const octave_idx_type ptype, FloatComplex *cout)
{
    dgtreal_ola_s(f, g, L, gl, W, a, M, bl,
                  static_cast<dgt_phasetype>(ptype),
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
    const octave_idx_type ptype = args(5).int_value();

    MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
    MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));

    const octave_idx_type L  = f.rows();
    const octave_idx_type W  = f.columns();
    const octave_idx_type gl = g.rows();
    const octave_idx_type N = L / a;

    const octave_idx_type M2 = M / 2 + 1;

    dim_vector dims_out(M2, N, W);
    dims_out.chop_trailing_singletons();

    MArray<LTFAT_COMPLEX> cout(dims_out);

    fwd_dgtreal_ola(f.data(), g.data(), L, gl, W, a, M, bl, ptype,
                    cout.fortran_vec());

    return octave_value(cout);
}
