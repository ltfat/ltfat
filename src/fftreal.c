#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

/**
* FFreal routines
*/

LTFAT_API void
LTFAT_NAME(fftreal)(LTFAT_REAL* f, ltfat_int L, ltfat_int W,
                    LTFAT_COMPLEX* cout)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(fftreal_init)(f, L, W, cout, FFTW_ESTIMATE);
    LTFAT_NAME(fftreal_execute)(p, f, cout);
    LTFAT_FFTW(destroy_plan)(p);
}

LTFAT_API void
LTFAT_NAME(fftreal_execute)(const LTFAT_FFTW(plan) p, LTFAT_REAL* f,
                            LTFAT_COMPLEX* cout)
{
    LTFAT_FFTW(execute_dft_r2c)(p, f, (LTFAT_FFTW(complex)*) cout);
}


/*
* IF anything else than FFTW_ESTIMATE is used for a flag, the planning overwrites input array !
*/
LTFAT_API LTFAT_FFTW(plan)
LTFAT_NAME(fftreal_init)(LTFAT_REAL* f, ltfat_int L, ltfat_int W,
                         LTFAT_COMPLEX* cout, unsigned flag)
{
    int Lint = (int) L;
    int L2int = (Lint / 2) + 1;
    LTFAT_FFTW(plan) p = LTFAT_FFTW(plan_many_dft_r2c)(1, &Lint,(int) W,
                         f, NULL, 1, Lint, (LTFAT_FFTW(complex)*) cout, NULL,
                         1, L2int, flag);

    return p;
}

LTFAT_API void
LTFAT_NAME(ifftreal)(LTFAT_COMPLEX* c, ltfat_int L, ltfat_int W,
                     LTFAT_REAL* f)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(ifftreal_init)(c, L, W, f, FFTW_ESTIMATE);
    LTFAT_NAME(ifftreal_execute)(p, c, L, W, f);
    LTFAT_FFTW(destroy_plan)(p);
}

LTFAT_API void
LTFAT_NAME(ifftreal_execute)(const LTFAT_FFTW(plan) p, LTFAT_COMPLEX* c,
                             ltfat_int L, ltfat_int W,
                             LTFAT_REAL* f)
{
    LTFAT_FFTW(execute_dft_c2r)(p, (LTFAT_FFTW(complex)*) c, f);

    LTFAT_REAL s  = (LTFAT_REAL) (1.0 / L);
    for (ltfat_int ii = 0; ii < L * W; ii++)
        f[ii] *= s;
}


/*
* IF anything else than FFTW_ESTIMATE is used for a flag, the planning overwrites input array !
*/
LTFAT_API LTFAT_FFTW(plan)
LTFAT_NAME(ifftreal_init)(LTFAT_COMPLEX* c, ltfat_int L, ltfat_int W,
                          LTFAT_REAL* f, unsigned flag)
{
    int Lint = (int) L;
    int L2int = (L / 2) + 1;
    LTFAT_FFTW(plan) p = LTFAT_FFTW(plan_many_dft_c2r)(1, &Lint,(int) W,
                         (LTFAT_FFTW(complex)*) c, NULL, 1, L2int,
                         f, NULL, 1, Lint, flag);

    return p;
}
