#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"


/**
* FFT filterbank routines
*/

LTFAT_EXTERN void
LTFAT_NAME(fftreal)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
                    LTFAT_COMPLEX *cout)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(fftreal_plan)(f, L, W, cout, FFTW_ESTIMATE);
    LTFAT_FFTW(execute_dft_r2c)(p,(LTFAT_REAL*)f,cout);
    LTFAT_FFTW(destroy_plan)(p);
}


/*
* IF anything else than FFTW_ESTIMATE is used for a flag, the planning overwrites input array !
*/
LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(fftreal_plan)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
                         LTFAT_COMPLEX *cout, unsigned flag)
{
    ltfatInt L2 = (L/2)+1;
    int ltfatInt = (int) L;
    LTFAT_FFTW(plan) p = LTFAT_FFTW(plan_many_dft_r2c)(1, &ltfatInt, W,
                         (LTFAT_REAL*)f, NULL,
                         1, L,
                         (LTFAT_COMPLEX*)cout, NULL,
                         1, L2,
                         flag);

    return p;
}

LTFAT_EXTERN void
LTFAT_NAME(ifftreal)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
                     LTFAT_REAL *f)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(ifftreal_plan)(c, L, W, f, FFTW_ESTIMATE);
    LTFAT_FFTW(execute_dft_c2r)(p,(LTFAT_COMPLEX*)c,f);

    LTFAT_REAL s  = (LTFAT_REAL) (1.0/L);
    for (ltfatInt ii=0; ii<L*W; ii++)
    {
        f[ii] *=s;
    }

    LTFAT_FFTW(destroy_plan)(p);
}


/*
* IF anything else than FFTW_ESTIMATE is used for a flag, the planning overwrites input array !
*/
LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(ifftreal_plan)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
                          LTFAT_REAL *f, unsigned flag)
{
    ltfatInt L2 = (L/2)+1;
    int ltfatInt = (int) L;
    LTFAT_FFTW(plan) p = LTFAT_FFTW(plan_many_dft_c2r)(1, &ltfatInt, W,
                         (LTFAT_COMPLEX*)c, NULL,
                         1, L2,
                         f, NULL,
                         1, L,
                         flag);

    return p;
}


