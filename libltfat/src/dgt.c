#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN
void LTFAT_NAME_COMPLEX(dgt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                                  const ltfatInt L, const ltfatInt W,
                                  const ltfatInt a, const ltfatInt M,
                                  const dgt_phasetype ptype, LTFAT_COMPLEX *cout)
{
    LTFAT_NAME(dgt_long_plan) plan =
        LTFAT_NAME(dgt_long_init)(f, g, L, W, a, M, cout, ptype, FFTW_ESTIMATE);

    LTFAT_NAME(dgt_long_execute)(plan);

    LTFAT_NAME(dgt_long_done)(plan);
}

LTFAT_EXTERN
void LTFAT_NAME(dgt_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                          const ltfatInt L, const ltfatInt W, const ltfatInt a,
                          const ltfatInt M, const dgt_phasetype ptype,
                          LTFAT_COMPLEX *cout)
{
    LTFAT_COMPLEX* gf = (LTFAT_COMPLEX*) ltfat_malloc(L * sizeof(LTFAT_COMPLEX));

    LTFAT_NAME_REAL(wfac)(g, L, 1, a, M, gf);

    LTFAT_NAME(dgt_fac_r)(f, (const LTFAT_COMPLEX*)gf, L, W, a, M, ptype, cout);

    ltfat_free(gf);
}


LTFAT_EXTERN
void LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                              const ltfatInt L, const ltfatInt W, const ltfatInt a,
                              const ltfatInt M, const dgt_phasetype ptype,
                              LTFAT_COMPLEX *cout)
{

    LTFAT_NAME(dgtreal_long_plan) plan =
        LTFAT_NAME(dgtreal_long_init)(f, g, L, W, a, M, cout, ptype, FFTW_ESTIMATE);

    LTFAT_NAME(dgtreal_long_execute)(plan);

    LTFAT_NAME(dgtreal_long_done)(plan);

}

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dgt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                           const ltfatInt L, const ltfatInt gl,
                           const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                           const dgt_phasetype ptype, LTFAT_COMPLEX *cout)
{

    LTFAT_NAME(dgt_fb_plan) plan =
        LTFAT_NAME(dgt_fb_init)(g, gl, a, M, ptype, FFTW_ESTIMATE);

    LTFAT_NAME(dgt_fb_execute)(plan, f, L, W, cout);

    LTFAT_NAME(dgt_fb_done)(plan);
}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                       const ltfatInt L, const ltfatInt gl,
                       const ltfatInt W, const ltfatInt a, const ltfatInt M,
                       const dgt_phasetype ptype, LTFAT_COMPLEX *cout)
{

    LTFAT_NAME(dgtreal_fb_plan) plan =
        LTFAT_NAME(dgtreal_fb_init)(g, gl, a, M, ptype, FFTW_ESTIMATE);

    LTFAT_NAME(dgtreal_fb_execute)(plan, f, L, W, cout);

    LTFAT_NAME(dgtreal_fb_done)(plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
                    const ltfatInt bl, const dgt_phasetype ptype,
                    LTFAT_COMPLEX *cout)
{
    LTFAT_NAME(dgt_ola_plan) plan =
        LTFAT_NAME(dgt_ola_init)(g, gl, W, a, M, bl, ptype, FFTW_ESTIMATE);

    LTFAT_NAME(dgt_ola_execute)(plan, f, L, cout);

    LTFAT_NAME(dgt_ola_done)(plan);

}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                        const ltfatInt L, const ltfatInt gl,
                        const ltfatInt W, const ltfatInt a, const ltfatInt M,
                        const ltfatInt bl, const dgt_phasetype ptype,
                        LTFAT_COMPLEX *cout)
{
    LTFAT_NAME(dgtreal_ola_plan) plan =
        LTFAT_NAME(dgtreal_ola_init)(g, gl, W, a, M, bl, ptype, FFTW_ESTIMATE);

    LTFAT_NAME(dgtreal_ola_execute)(plan, f, L, cout);

    LTFAT_NAME(dgtreal_ola_done)(plan);

}


