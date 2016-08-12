#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME(wfac_plan)
{
    ltfatInt b;
    ltfatInt c;
    ltfatInt p;
    ltfatInt q;
    ltfatInt d;
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    LTFAT_REAL scaling;
    LTFAT_REAL* sbuf;
    LTFAT_FFTW(plan) p_before;
};

LTFAT_EXTERN int
LTFAT_NAME(wfac_init)(const ltfatInt L, const ltfatInt a, const ltfatInt M,
                      unsigned flags, LTFAT_NAME(wfac_plan)** pout)
{
    LTFAT_NAME(wfac_plan)* plan = NULL;
    ltfatInt minL, h_a, h_m;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(pout);
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a (passed %d) must be positive.", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M (passed %d) must be positive.", M);

    minL = ltfat_lcm(a, M);
    CHECK(LTFATERR_BADARG,
          L > 0 && !(L % minL),
          "L (passed %d) must be positive and divisible by lcm(a,M)=%d.",
          L, minL);

    CHECKMEM(plan = (LTFAT_NAME(wfac_plan)*) ltfat_calloc(1, sizeof * plan));

    plan->b = L / M;
    plan->c = ltfat_gcd(a, M, &h_a, &h_m);
    plan->p = a / plan->c;
    plan->q = M / plan->c;
    plan->d = plan->b / plan->p;
    plan->scaling = (LTFAT_REAL) ( sqrt((double)M) );
    plan->a = a; plan->M = M; plan->L = L;

    CHECKMEM(plan->sbuf = LTFAT_NAME_REAL(malloc)(2 * plan->d));

    /* Create plan. In-place. */
    plan->p_before = LTFAT_FFTW(plan_dft_1d)((int)plan->d,
                     (LTFAT_FFTW(complex)*) plan->sbuf,
                     (LTFAT_FFTW(complex)*) plan->sbuf,
                     FFTW_FORWARD, flags);

    CHECKINIT(plan->p_before, "FFTW plan creation failed.");

    *pout = plan;
    return status;
error:
    if (plan)
    {
        if (plan->p_before) LTFAT_FFTW(destroy_plan)(plan->p_before);
        ltfat_free(plan->sbuf);
        ltfat_free(plan);
    }
    *pout = NULL;
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(wfac_execute)(LTFAT_NAME(wfac_plan)* plan, const LTFAT_TYPE* g,
                         const ltfatInt R, LTFAT_COMPLEX* gf)
{
    ltfatInt rem, negrem, c, p, q, d, a, M, L, ld3;
    LTFAT_REAL* sbuf, *gfp;
    LTFAT_REAL scaling;
    LTFAT_FFTW(plan) p_before;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(g); CHECKNULL(gf);
    CHECK(LTFATERR_NOTPOSARG, R > 0, "R (passed %d) must be positive.", R);

    sbuf = plan->sbuf;
    p_before = plan->p_before;

    /* const ltfatInt b = plan->b; */
    c = plan->c;
    p = plan->p;
    q = plan->q;
    d = plan->d;
    scaling = plan->scaling;
    a = plan->a;
    M = plan->M;
    L = plan->L;

    ld3 = c * p * q * R;
    gfp = (LTFAT_REAL*)gf;

    for (ltfatInt r = 0; r < c; r++)
    {
        for (ltfatInt w = 0; w < R; w++)
        {
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt k = 0; k < p; k++)
                {
                    negrem = ltfat_positiverem(k * M - l * a, L);
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        rem = (negrem + s * p * M) % L;
#ifdef LTFAT_COMPLEXTYPE
                        LTFAT_COMPLEX gval = scaling * g[r + rem + L * w];
                        sbuf[2 * s]   = ltfat_real(gval);
                        sbuf[2 * s + 1] = ltfat_imag(gval);
#else
                        sbuf[2 * s]   = scaling * g[r + rem + L * w];
                        sbuf[2 * s + 1] = 0.0;
#endif
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s = 0; s < 2 * d; s += 2)
                    {
                        gfp[s * ld3]  = sbuf[s];
                        gfp[s * ld3 + 1] = sbuf[s + 1];
                    }
                    gfp += 2;
                }
            }
        }
    }

error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(wfac_done)(LTFAT_NAME(wfac_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(pout);
    CHECKNULL(*pout);

    LTFAT_FFTW(destroy_plan)((*pout)->p_before);
    ltfat_free((*pout)->sbuf);
    ltfat_free(*pout);
    *pout = NULL;
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(wfac)(const LTFAT_TYPE* g, const ltfatInt L, const ltfatInt R,
                 const ltfatInt a, const ltfatInt M,
                 LTFAT_COMPLEX* gf)
{
    LTFAT_NAME(wfac_plan)* p = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKSTATUS( LTFAT_NAME(wfac_init)( L, a, M, FFTW_ESTIMATE, &p),
                 "Init failed");

    CHECKSTATUS( LTFAT_NAME(wfac_execute)(p, g, R, gf), "Execute failed");

error:
    if (p) LTFAT_NAME(wfac_done)(&p);
    return status;
}


/* LTFAT_EXTERN int */
/* LTFAT_NAME(wfac)(const LTFAT_TYPE* g, const ltfatInt L, const ltfatInt R, */
/*                  const ltfatInt a, const ltfatInt M, */
/*                  LTFAT_COMPLEX* gf) */
/* { */
/*     ltfatInt h_a, h_m, s; */
/*  */
/*     LTFAT_REAL* sbuf, *gfp; */
/*  */
/*     ltfatInt rem, negrem; */
/*  */
/*     LTFAT_FFTW(plan) p_before; */
/*  */
/*     const ltfatInt b = L / M; */
/*     const ltfatInt c = ltfat_gcd(a, M, &h_a, &h_m); */
/*     const ltfatInt p = a / c; */
/*     const ltfatInt q = M / c; */
/*     const ltfatInt d = b / p; */
/*  */
/*     const double sqrtM = sqrt(M); */
/*  */
/*     sbuf = ltfat_malloc(2 * d * sizeof * sbuf); */
/*  */
/*     #<{(| Create plan. In-place. |)}># */
/*     p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, */
/*                                        (LTFAT_COMPLEX*)sbuf, */
/*                                        FFTW_FORWARD, FFTW_MEASURE); */
/*  */
/*     const ltfatInt ld3 = c * p * q * R; */
/*     gfp = (LTFAT_REAL*)gf; */
/*     for (ltfatInt r = 0; r < c; r++) */
/*     { */
/*         for (ltfatInt w = 0; w < R; w++) */
/*         { */
/*             for (ltfatInt l = 0; l < q; l++) */
/*             { */
/*                 for (ltfatInt k = 0; k < p; k++) */
/*                 { */
/*                     negrem = ltfat_positiverem(k * M - l * a, L); */
/*                     for (s = 0; s < d; s++) */
/*                     { */
/*                         rem = (negrem + s * p * M) % L; */
/* #ifdef LTFAT_COMPLEXTYPE */
/*                         LTFAT_COMPLEX gval = sqrtM * g[r + rem + L * w]; */
/*                         sbuf[2 * s]   = creal(gval); */
/*                         sbuf[2 * s + 1] = cimag(gval); */
/* #else */
/*                         sbuf[2 * s]   = sqrtM * g[r + rem + L * w]; */
/*                         sbuf[2 * s + 1] = 0.0; */
/* #endif */
/*                     } */
/*  */
/*                     LTFAT_FFTW(execute)(p_before); */
/*  */
/*                     for (s = 0; s < 2 * d; s += 2) */
/*                     { */
/*                         gfp[s * ld3]  = sbuf[s]; */
/*                         gfp[s * ld3 + 1] = sbuf[s + 1]; */
/*                     } */
/*                     gfp += 2; */
/*                 } */
/*             } */
/*         } */
/*     } */
/*  */
/*     ltfat_free(sbuf); */
/*     LTFAT_FFTW(destroy_plan)(p_before); */
/* } */
