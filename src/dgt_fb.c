#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME(dgt_fb_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt gl;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_small;
    LTFAT_COMPLEX* sbuf;
    LTFAT_COMPLEX* fw;
    LTFAT_TYPE* gw;
};


LTFAT_EXTERN int
LTFAT_NAME(dgt_fb)(const LTFAT_TYPE* f, const LTFAT_TYPE* g,
                   const ltfatInt L, const ltfatInt gl,
                   const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                   const dgt_phasetype ptype, LTFAT_COMPLEX* cout)
{

    LTFAT_NAME(dgt_fb_plan)* plan = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        LTFAT_NAME(dgt_fb_init)(g, gl, a, M, ptype, FFTW_ESTIMATE, &plan),
        "Init failed");

    CHECKSTATUS(
        LTFAT_NAME(dgt_fb_execute)(plan, f, L, W, cout),
        "Execute failed");

error:
    if (plan) LTFAT_NAME(dgt_fb_done)(&plan);
    return status;
}

/* The following macro adds the coefficients together performing the
 * last part of the Poisson summation, executes the FFT on the summed
 * coefficients, and places the coefficients in the output array.
 *
 * The first summation is done in that peculiar way to obtain the
 * correct phase for a frequency invariant Gabor transform. Summing
 * them directly would lead to a time invariant (phase-locked) Gabor
 * transform.
 *
 * The macro is called in three different places in the dgt_fb function.
 */
#define THE_SUM { \
LTFAT_NAME_COMPLEX(fold_array)(fw,gl,plan.ptype?-glh:n*a-glh,M,sbuf); \
LTFAT_FFTW(execute)(plan.p_small); \
memcpy(cout + (n*M + w*M*N),sbuf,M*sizeof*cout); \
}


LTFAT_EXTERN int //LTFAT_NAME(dgt_fb_plan)
LTFAT_NAME(dgt_fb_init)(const LTFAT_TYPE* g,
                        const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, unsigned flags, LTFAT_NAME(dgt_fb_plan)** p)
{
    LTFAT_NAME(dgt_fb_plan)* plan = NULL;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(g);
    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    CHECKMEM(plan = ltfat_calloc(1, sizeof * plan));

    plan->a = a;
    plan->M = M;
    plan->gl = gl;
    plan->ptype = ptype;

    CHECKMEM(plan->gw  = ltfat_malloc(plan->gl * sizeof * plan->gw));
    CHECKMEM(plan->fw  = ltfat_calloc(plan->gl, sizeof * plan->fw));
    CHECKMEM(plan->sbuf = ltfat_malloc(M * sizeof * plan->sbuf));

    plan->p_small = LTFAT_FFTW(plan_dft_1d)(M, (LTFAT_COMPLEX*)plan->sbuf,
                                            (LTFAT_COMPLEX*)plan->sbuf,
                                            FFTW_FORWARD, flags);

    CHECKINIT(plan->p_small, "FFTW plan creation failed.");

    LTFAT_NAME(fftshift)(g, gl, plan->gw);
    LTFAT_NAME(conjugate_array)(plan->gw, gl, plan->gw);

    // Assign the "return" value
    *p = plan;
    return status;
error:
    if (plan)
    {
        if (plan->p_small) LTFAT_FFTW(destroy_plan)(plan->p_small);
        LTFAT_SAFEFREEALL(plan->gw, plan->fw, plan->sbuf);
        ltfat_free(plan);
    }
    *p = NULL;
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(dgt_fb_done)(LTFAT_NAME(dgt_fb_plan)** plan)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(*plan);
    LTFAT_NAME(dgt_fb_plan)* pp = *plan;

    LTFAT_SAFEFREEALL(pp->sbuf, pp->gw, pp->fw);
    LTFAT_FFTW(destroy_plan)(pp->p_small);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(dgt_fb_execute)(const LTFAT_NAME(dgt_fb_plan)* p,
                           const LTFAT_TYPE* f,
                           const ltfatInt L, const ltfatInt W,  LTFAT_COMPLEX* cout)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(f); CHECKNULL(cout);
    CHECK(LTFATERR_BADARG, L >= p->gl && !(L % p->a) ,
          "L (passed %d) must be positive and divisible by a (passed %d).", L, p->a);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");

    /*  --------- initial declarations -------------- */
    LTFAT_NAME(dgt_fb_plan) plan = *p;

    const ltfatInt a = plan.a;
    const ltfatInt M = plan.M;
    const ltfatInt N = L / a;

    const ltfatInt gl = plan.gl;
    LTFAT_COMPLEX* sbuf = plan.sbuf;
    LTFAT_COMPLEX* fw = plan.fw;

    /* This is a floor operation. */
    const ltfatInt glh = plan.gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));

    LTFAT_TYPE* fbd;

    /*  ---------- main body ----------- */

    /*----- Handle the first boundary using periodic boundary conditions.*/
    for (ltfatInt n = 0; n < glh_d_a; n++)
    {
        for (ltfatInt w = 0; w < W; w++)
        {

            fbd = (LTFAT_TYPE*)f + (L - (glh - n * a) + L * w);
            for (ltfatInt l = 0; l < glh - n * a; l++)
                fw[l] = fbd[l] * plan.gw[l];

            fbd = (LTFAT_TYPE*)f -  (glh - n * a) +  L * w;
            for (ltfatInt l = glh - n * a; l < gl; l++)
                fw[l] = fbd[l] * plan.gw[l];

            THE_SUM

        }
    }

    /* ----- Handle the middle case. --------------------- */
    for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = (LTFAT_TYPE*)f + (n * a - glh + L * w);
            for (ltfatInt l = 0; l < gl; l++)
                fw[l] = fbd[l] * plan.gw[l];

            THE_SUM
        }

    }

    /* Handle the last boundary using periodic boundary conditions. */
    for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = (LTFAT_TYPE*)f + (n * a - glh + L * w);
            for (ltfatInt l = 0; l < L - n * a + glh; l++)
                fw[l] = fbd[l] * plan.gw[l];

            fbd = (LTFAT_TYPE*)f - (L - n * a + glh) +  L * w;
            for (ltfatInt l = L - n * a + glh; l < gl; l++)
                fw[l] = fbd[l] * plan.gw[l];

            THE_SUM
        }
    }

error:
    return status;
}

#undef THE_SUM
