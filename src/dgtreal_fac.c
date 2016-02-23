#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN LTFAT_NAME(dgtreal_long_plan)
LTFAT_NAME(dgtreal_long_init)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                              const ltfatInt L, const ltfatInt W,
                              const ltfatInt a, const ltfatInt M,
                              LTFAT_COMPLEX *cout, const dgt_phasetype ptype,
                              unsigned flags)
{

    LTFAT_NAME(dgtreal_long_plan) plan;
    ltfatInt h_m;

    plan.a = a;
    plan.M = M;
    plan.L = L;
    plan.W = W;
    plan.ptype = ptype;
    const ltfatInt N = L / a;



    plan.c = gcd(a, M, &plan.h_a, &h_m);
    const ltfatInt b = L / M;
    const ltfatInt p = a / plan.c;
    const ltfatInt q = M / plan.c;
    const ltfatInt d = b / p;
    plan.h_a = -plan.h_a;

    const ltfatInt M2 = M / 2 + 1;
    const ltfatInt d2 = d / 2 + 1;

    plan.sbuf = ltfat_malloc( d * sizeof(LTFAT_REAL));
    plan.cbuf = ltfat_malloc(d2 * sizeof(LTFAT_COMPLEX));
    plan.cout = cout;
    plan.f    = f;

    plan.ff = ltfat_malloc(2 * d2 * p * q * W * sizeof(LTFAT_REAL));
    plan.cf = ltfat_malloc(2 * d2 * q * q * W * sizeof(LTFAT_REAL));

    const ltfatInt wfs = wfacreal_size(L, a, M);

    plan.gf   = (LTFAT_COMPLEX*)ltfat_malloc(wfs * sizeof(LTFAT_COMPLEX));

    plan.cwork = (LTFAT_REAL*)ltfat_malloc(M * N * W * sizeof(LTFAT_REAL));

    /* Get factorization of window */
    LTFAT_NAME(wfacreal)(g, L, 1, a, M, plan.gf);

    /* Create plans. In-place. */
    // Downcast to int
    int Mint = (int) plan.M;

    plan.p_veryend =
        LTFAT_FFTW(plan_many_dft_r2c)(1, &Mint, N * W,
                                      plan.cwork, NULL,
                                      1, M,
                                      cout, NULL,
                                      1, M2,
                                      flags);

    plan.p_before =
        LTFAT_FFTW(plan_dft_r2c_1d)(d, plan.sbuf, plan.cbuf, flags);

    plan.p_after  =
        LTFAT_FFTW(plan_dft_c2r_1d)(d, plan.cbuf, plan.sbuf, flags);

    return plan;
}




LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long_execute)(const LTFAT_NAME(dgtreal_long_plan) plan)
{

    LTFAT_NAME(dgtreal_walnut_plan)(plan);

    if (plan.ptype)
    {
        LTFAT_NAME_REAL(dgtphaselockhelper)(plan.cwork, plan.L, plan.W, plan.a,
                                            plan.M, plan.cwork);
        /*ltfatInt N = plan.L / plan.a;
        ltfatInt W = plan.W;
        ltfatInt M = plan.M;

        for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt n = 0; n < N; n++)
            {
                LTFAT_REAL* cworktmp = plan.cwork + w * N * M + n * M;
                LTFAT_NAME_REAL(circshift)(cworktmp, cworktmp, M, -plan.a * n);
            }

        }
        */
    }

    /* FFT to modulate the coefficients. */
    LTFAT_FFTW(execute)(plan.p_veryend);

}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan) plan)
{
    LTFAT_FFTW(destroy_plan)(plan.p_veryend);
    LTFAT_FFTW(destroy_plan)(plan.p_before);
    LTFAT_FFTW(destroy_plan)(plan.p_after);
    LTFAT_SAFEFREEALL(plan.sbuf, plan.cbuf, plan.cwork, plan.gf, plan.ff, plan.cf);
}
