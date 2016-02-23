#include "ltfat.h"
#include "ltfat_types.h"


/**
* FFreal routines
*/

LTFAT_EXTERN void
LTFAT_NAME(fftreal)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
                    LTFAT_COMPLEX *cout)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(fftreal_init)(f, L, W, cout, FFTW_ESTIMATE);
    LTFAT_NAME(fftreal_execute)(p,f,cout);
    LTFAT_FFTW(destroy_plan)(p);
}

LTFAT_EXTERN void
LTFAT_NAME(fftreal_execute)(const LTFAT_FFTW(plan) p, LTFAT_REAL *f,
                            LTFAT_COMPLEX *cout)
{
    LTFAT_FFTW(execute_dft_r2c)(p,f,cout);
}


/*
* IF anything else than FFTW_ESTIMATE is used for a flag, the planning overwrites input array !
*/
LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(fftreal_init)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
                         LTFAT_COMPLEX *cout, unsigned flag)
{
    ltfatInt L2 = (L/2)+1;
    int ltfatInt = (int) L;
    LTFAT_FFTW(plan) p = LTFAT_FFTW(plan_many_dft_r2c)(1, &ltfatInt, W,
                         f, NULL,
                         1, L,
                         cout, NULL,
                         1, L2,
                         flag);

    return p;
}

LTFAT_EXTERN void
LTFAT_NAME(ifftreal)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
                     LTFAT_REAL *f)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(ifftreal_init)(c, L, W, f, FFTW_ESTIMATE);
    LTFAT_NAME(ifftreal_execute)(p,c,L,W,f);
    LTFAT_FFTW(destroy_plan)(p);
}

LTFAT_EXTERN void
LTFAT_NAME(ifftreal_execute)(const LTFAT_FFTW(plan) p, LTFAT_COMPLEX *c,
                             const ltfatInt L, const ltfatInt W,
                             LTFAT_REAL *f)
{
    LTFAT_FFTW(execute_dft_c2r)(p,c,f);

    LTFAT_REAL s  = (LTFAT_REAL) (1.0/L);
    for (ltfatInt ii=0; ii<L*W; ii++)
    {
        f[ii] *=s;
    }
}


/*
* IF anything else than FFTW_ESTIMATE is used for a flag, the planning overwrites input array !
*/
LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(ifftreal_init)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
                          LTFAT_REAL *f, unsigned flag)
{
    ltfatInt L2 = (L/2)+1;
    int ltfatInt = (int) L;
    LTFAT_FFTW(plan) p = LTFAT_FFTW(plan_many_dft_c2r)(1, &ltfatInt, W,
                         c, NULL,
                         1, L2,
                         f, NULL,
                         1, L,
                         flag);

    return p;
}


