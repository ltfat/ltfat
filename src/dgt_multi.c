#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

/* struct LTFAT_NAME(dgt_multi_plan) */
/* { */
/*     ltfatInt a; */
/*     ltfatInt M; */
/*     ltfatInt L; */
/*     ltfatInt Lg; */
/*     ltfatInt W; */
/*     ltfatInt lt1; */
/*     ltfatInt lt2; */
/*  */
/*     LTFAT_COMPLEX *f; */
/*     LTFAT_COMPLEX *c_scratch; */
/*     LTFAT_COMPLEX *cout; */
/*  */
/*     LTFAT_COMPLEX *mwin; */
/*     LTFAT_COMPLEX *c_rect; */
/*  */
/*     LTFAT_COMPLEX *mod; */
/*  */
/*     LTFAT_NAME(dgt_long_plan)** rect_plan_array; */
/* }; */

LTFAT_EXTERN void
LTFAT_NAME(nonsepwin2multi)(const LTFAT_COMPLEX *g,
                            const ltfatInt L, const ltfatInt Lg, const ltfatInt a, const ltfatInt M,
                            const ltfatInt lt1, const ltfatInt lt2,
                            LTFAT_COMPLEX *mwin)
{
    const ltfatInt b = L / M;

    const LTFAT_REAL scal = 2 * M_PI / L;

    LTFAT_COMPLEX *gwork = (LTFAT_COMPLEX *)ltfat_malloc(L * sizeof(LTFAT_COMPLEX));
    LTFAT_NAME_COMPLEX(fir2long)((const LTFAT_COMPLEX*)g, Lg, L, (LTFAT_COMPLEX*)gwork);

    for (ltfatInt w = 0; w < lt2; w++)
    {
        const ltfatInt wavenum = ((w * lt1) % lt2) * b / lt2;
        for (ltfatInt l = 0; l < L; l++)
        {
            mwin[l + w * L] = cexp(I * scal * l * wavenum) * gwork[positiverem(l - w * a, L)];
        }
    }

    ltfat_free(gwork);
}


LTFAT_EXTERN LTFAT_NAME(dgt_multi_plan)
LTFAT_NAME(dgt_multi_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                           const ltfatInt L, const ltfatInt Lg,
                           const ltfatInt W, const ltfatInt a, const ltfatInt M,
                           const ltfatInt lt1, const ltfatInt lt2,
                           LTFAT_COMPLEX *cout, unsigned flags)
{

    LTFAT_NAME(dgt_multi_plan) plan;

    plan.a = a;
    plan.M = M;
    plan.L = L;
    plan.Lg = Lg;
    plan.W = W;

    plan.lt1 = lt1;
    plan.lt2 = lt2;

    plan.f     = (LTFAT_COMPLEX *)f;
    plan.cout  = cout;

    const ltfatInt N   = L / a;
    const ltfatInt Ns  = N / lt2;

    plan.mwin = ltfat_malloc(L * lt2 * sizeof(LTFAT_COMPLEX));
    plan.c_scratch = ltfat_malloc(M * Ns * W * sizeof(LTFAT_COMPLEX));


    LTFAT_NAME(nonsepwin2multi)(g, L, Lg, a, M, lt1, lt2, plan.mwin);

    plan.rect_plan_array = ltfat_malloc(lt2 * sizeof(LTFAT_NAME_COMPLEX(dgt_long_plan)*));

    for (ltfatInt win = 0; win < lt2; win++)
    {
            LTFAT_NAME_COMPLEX(dgt_long_init)(plan.f,
                                      plan.mwin + L * win,
                                      L, W, a * lt2, M, 
                                      plan.c_scratch, 0, flags,
                                      &plan.rect_plan_array[win]);
    }

    plan.mod = (LTFAT_COMPLEX*) ltfat_malloc(N * sizeof(LTFAT_COMPLEX));

    for (ltfatInt win = 0; win < plan.lt2; win++)
    {
        for (ltfatInt n = 0; n < Ns; n++)
        {
            plan.mod[win + n * lt2] = cexp(-2 * M_PI * I * (a * n * ((win * lt1) % lt2)) / M);
        }
    }

    return plan;
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_multi_execute)(const LTFAT_NAME(dgt_multi_plan) plan)
{
    const ltfatInt N   = plan.L / plan.a;
    const ltfatInt Ns  = N / plan.lt2;

    const ltfatInt M = plan.M;
    const ltfatInt W = plan.W;
    const ltfatInt lt2 = plan.lt2;

    for (ltfatInt win = 0; win < plan.lt2; win++)
    {
        LTFAT_NAME_COMPLEX(dgt_long_execute)(plan.rect_plan_array[win]);
        for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt n = 0; n < Ns; n++)
            {
                for (ltfatInt m = 0; m < M; m++)
                {
                    plan.cout[m + win * M + n * lt2 * M + w * M * N] = plan.mod[win + n * lt2] * plan.c_scratch[m + n * M + w * M * Ns];
                }
            }
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_multi_done)(LTFAT_NAME(dgt_multi_plan) plan)
{
    for (ltfatInt ii = 0; ii < plan.lt2; ii++)
    {
        LTFAT_NAME_COMPLEX(dgt_long_done)(&plan.rect_plan_array[ii]);
    }
    LTFAT_SAFEFREEALL(plan.mod, plan.rect_plan_array, plan.c_scratch, plan.mwin);
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_multi)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                      const ltfatInt L, const ltfatInt Lg, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                      const ltfatInt lt1, const ltfatInt lt2,
                      LTFAT_COMPLEX *cout)
{

    LTFAT_NAME(dgt_multi_plan) plan = LTFAT_NAME(dgt_multi_init)(
                                          f, g, L, Lg, W, a, M, lt1, lt2, cout, FFTW_ESTIMATE);

    LTFAT_NAME(dgt_multi_execute)(plan);

    LTFAT_NAME(dgt_multi_done)(plan);

}
