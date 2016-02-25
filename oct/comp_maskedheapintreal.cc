#define TYPEDEPARGS 0
#define MATCHEDARGS 1, 2, 8
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_maskedheapintreal
#define OCTFILEHELP "Computes masked heapint.\n\
Usage: c = comp_maskedheapintreal(s, itime, ifreq, mask, a, M, tol, do_timeinv, usephase);\n\
Yeah."


#include "ltfat_oct_template_helper.h"

static inline void
fwd_maskedheapintreal(const double *s, const double *tgrad, const double *fgrad,
                      const int* mask, const octave_idx_type a, const octave_idx_type M,
                      const octave_idx_type L, const octave_idx_type W,
                      const double tol, dgt_phasetype phasetype,  double *phase)
{
    if (phasetype == 2)
        maskedheapintreal_relgrad_d( s, tgrad, fgrad, mask, a, M, L, W, tol, phasetype, phase);
    else
        maskedheapintreal_d( s, tgrad, fgrad, mask, a, M, L, W, tol, phase);
}

static inline void
fwd_maskedheapintreal(const float *s, const float *tgrad, const float *fgrad,
                      const int* mask, const octave_idx_type a, const octave_idx_type M,
                      const octave_idx_type L, const octave_idx_type W,
                      const float tol, dgt_phasetype phasetype,  float *phase)
{
    if (phasetype == 2)
        maskedheapintreal_relgrad_s( s, tgrad, fgrad, mask, a, M, L, W, tol, phasetype, phase);
    else
        maskedheapintreal_s( s, tgrad, fgrad, mask, a, M, L, W, tol, phase);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    MArray<LTFAT_REAL> s = ltfatOctArray<LTFAT_REAL>(args(0));
    MArray<LTFAT_REAL> tgrad = ltfatOctArray<LTFAT_REAL>(args(1));
    MArray<LTFAT_REAL> fgrad = ltfatOctArray<LTFAT_REAL>(args(2));
    const int32NDArray mask  = args(3).int32_array_value();
    const octave_idx_type a  = args(4).int_value();
    const octave_idx_type M  = args(5).int_value();
    const double tol  = args(6).double_value();
    const octave_idx_type phasetype  = args(7).int_value();
    MArray<LTFAT_REAL> usephase = ltfatOctArray<LTFAT_REAL>(args(8));

    const octave_idx_type M2 = args(0).rows();
    const octave_idx_type N = args(0).columns();
    const octave_idx_type L = N * a;
    const octave_idx_type W = s.nelem() / (M * N);

    MArray<LTFAT_REAL> phase(dim_vector(M2, N, W));

    memcpy(phase.fortran_vec(), usephase.data(), M2 * N * W * sizeof(LTFAT_REAL));

    fwd_maskedheapintreal(s.data(), tgrad.data(), fgrad.data(),
                          reinterpret_cast<const int*>(mask.data()),
                          a, M, L, W, static_cast<LTFAT_REAL>(tol),
                          static_cast<dgt_phasetype>(phasetype),
                          phase.fortran_vec());

    return octave_value(phase);
}
