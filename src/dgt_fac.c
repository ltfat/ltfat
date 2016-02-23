#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN LTFAT_NAME(dgt_long_plan)
LTFAT_NAME(dgt_long_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *cout,
                          const dgt_phasetype ptype, unsigned flags)
{

    LTFAT_NAME(dgt_long_plan) plan;
    ltfatInt h_m;

    plan.a = a;
    plan.M = M;
    plan.L = L;
    plan.W = W;
    plan.ptype = ptype;
    const ltfatInt N = L / a;
    const ltfatInt b = L / M;

    plan.c = gcd(a, M, &plan.h_a, &h_m);
    const ltfatInt p = a / plan.c;
    const ltfatInt q = M / plan.c;
    const ltfatInt d = b / p;
    plan.h_a = -plan.h_a;

    plan.sbuf = ltfat_malloc(2 * d * sizeof(LTFAT_REAL));
    plan.cout = cout;
    plan.f    = f;

    plan.gf   = ltfat_malloc(L * sizeof(LTFAT_COMPLEX));

    plan.ff = ltfat_malloc(2 * d * p * q * W * sizeof(LTFAT_REAL));
    plan.cf = ltfat_malloc(2 * d * q * q * W * sizeof(LTFAT_REAL));

    /* Get factorization of window */
    LTFAT_NAME_COMPLEX(wfac)(g, L, 1, a, M, plan.gf);

    // Downcasting to ints
    int Mint = (int) M;

    /* Create plans. In-place. */
    plan.p_veryend =
        LTFAT_FFTW(plan_many_dft)(1, &Mint, N * W,
                                  plan.cout, NULL,
                                  1, Mint,
                                  plan.cout, NULL,
                                  1, Mint,
                                  FFTW_FORWARD, flags);

    plan.p_before =
        LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)plan.sbuf,
                                (LTFAT_COMPLEX*)plan.sbuf,
                                FFTW_FORWARD, flags);

    plan.p_after  =
        LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)plan.sbuf,
                                (LTFAT_COMPLEX*)plan.sbuf,
                                FFTW_BACKWARD, flags);

    return plan;
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_long_execute)(const LTFAT_NAME(dgt_long_plan) plan)
{

    LTFAT_NAME(dgt_walnut_plan)(plan);

    if (plan.ptype)
    {

        LTFAT_NAME_COMPLEX(dgtphaselockhelper)(plan.cout, plan.L, plan.W, plan.a,
                                               plan.M, plan.cout);
        /*ltfatInt N = plan.L / plan.a;
        ltfatInt W = plan.W;

        for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt n = 0; n < N; n++)
            {
                LTFAT_COMPLEX* couttmp = plan.cout + w * N * plan.M + n * plan.M;
                LTFAT_NAME_COMPLEX(circshift)(couttmp, couttmp, plan.M, -plan.a * n);
            }

        }
        */
    }


    /* FFT to modulate the coefficients. */
    LTFAT_FFTW(execute)(plan.p_veryend);

}


LTFAT_EXTERN void
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan) plan)
{
    LTFAT_FFTW(destroy_plan)(plan.p_veryend);
    LTFAT_FFTW(destroy_plan)(plan.p_before);
    LTFAT_FFTW(destroy_plan)(plan.p_after);
    LTFAT_SAFEFREEALL(plan.sbuf, plan.gf, plan.ff, plan.cf);
}





/* -------------- Real valued signal ------------------------ */

LTFAT_EXTERN
void LTFAT_NAME(dgt_fac_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
                           const ltfatInt L, const ltfatInt W,
                           const ltfatInt a, const ltfatInt M,
                           const dgt_phasetype ptype, LTFAT_COMPLEX *cout)
{

    const ltfatInt N = L / a;
    // Downcasting to int
    int Mint = (int) M;

    /* Create plan. In-place. */
    LTFAT_FFTW(plan) p_veryend =
        LTFAT_FFTW(plan_many_dft)(1, &Mint, N * W,
                                  cout, NULL,
                                  1, Mint,
                                  cout, NULL,
                                  1, Mint,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_NAME(dgt_walnut_r)(f, gf, L, W, a, M, cout);


    if (ptype)
    {
        LTFAT_NAME_COMPLEX(dgtphaselockhelper)(cout, L, W, a, M, cout);
    }


    /* FFT to modulate the coefficients. */
    LTFAT_FFTW(execute)(p_veryend);

    LTFAT_FFTW(destroy_plan)(p_veryend);

}
