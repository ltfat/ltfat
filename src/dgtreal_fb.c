#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME(dgtreal_fb_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt gl;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_small;
    LTFAT_REAL*    sbuf;
    LTFAT_COMPLEX* cbuf;
    LTFAT_REAL* fw;
    LTFAT_REAL* gw;
    LTFAT_COMPLEX* cout;
};

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL* f, const LTFAT_REAL* g,
                       const ltfatInt L, const ltfatInt gl,
                       const ltfatInt W, const ltfatInt a, const ltfatInt M,
                       const dgt_phasetype ptype, LTFAT_COMPLEX* cout)
{
    LTFAT_NAME(dgtreal_fb_plan)* plan = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        LTFAT_NAME(dgtreal_fb_init)(g, gl, a, M, ptype, FFTW_ESTIMATE, &plan),
        "Init failed");

    CHECKSTATUS(
        LTFAT_NAME(dgtreal_fb_execute)(plan, f, L, W, cout),
        "Execute failed");

error:
    if (plan) LTFAT_NAME(dgtreal_fb_done)(&plan);
    return status;
}

/* #define THE_SUM_REAL { \ */
/* ltfatInt premarg = plan.ptype?-glh:n*a-glh; \ */
/* for (ltfatInt m=0;m<M;m++) \ */
/* { \ */
/*    rem = ltfat_positiverem(m+(premarg), M); \ */
/*    sbuf[rem]=0.0; \ */
/*    fbd=fw+m; \ */
/*    for (ltfatInt k=0;k<gl/M;k++) \ */
/*    { \ */
/*      sbuf[rem]+=(*fbd);         \ */
/*       fbd+=M; \ */
/*    } \ */
/* } \ */
/* \ */
/*  LTFAT_FFTW(execute)(plan.p_small);         \ */
/* \ */
/* coefsum=(LTFAT_REAL*)cout+2*(n*M2+w*M2*N); \ */
/* for (ltfatInt m=0;m<M2;m++) \ */
/* { \ */
/*    coefsum[2*m]   = CH(creal)(cbuf[m]); \ */
/*    coefsum[2*m+1] = CH(cimag)(cbuf[m]); \ */
/* }} */

#define THE_SUM_REAL { \
LTFAT_NAME(fold_array)(fw,gl,plan->ptype?-glh:n*a-glh,M,sbuf); \
LTFAT_FFTW(execute)(plan->p_small); \
memcpy(cout+(n*M2+w*M2*N),cbuf,M2*sizeof*cbuf); \
}


LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL* g,
                            const ltfatInt gl, const ltfatInt a,
                            const ltfatInt M, const dgt_phasetype ptype,
                            unsigned flags, LTFAT_NAME(dgtreal_fb_plan)** pout)
{
    LTFAT_NAME(dgtreal_fb_plan)* plan = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g);
    CHECKNULL(pout);
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    CHECKMEM(plan = ltfat_calloc(1, sizeof * plan));

    plan->a = a;
    plan->M = M;
    const ltfatInt M2 = M / 2 + 1;
    plan->gl = gl;
    plan->ptype = ptype;

    CHECKMEM( plan->gw   = ltfat_malloc(gl * sizeof * plan->gw));
    CHECKMEM( plan->fw   = ltfat_malloc(gl * sizeof * plan->fw));
    CHECKMEM( plan->sbuf = ltfat_malloc(M  * sizeof * plan->sbuf));
    CHECKMEM( plan->cbuf = ltfat_malloc(M2 * sizeof * plan->cbuf));

    plan->p_small = LTFAT_FFTW(plan_dft_r2c_1d)(M, plan->sbuf, plan->cbuf, flags);

    CHECKINIT(plan->p_small, "FFTW plan failed.");

    /* This is a floor operation. */
    /* const ltfatInt glh = gl / 2; */

    LTFAT_NAME(fftshift)(g, gl, plan->gw);
    /* for (ltfatInt l = 0; l < glh; l++) */
    /* { */
    /*     plan.gw[l] = g[l + (gl - glh)]; */
    /* } */
    /* for (ltfatInt l = glh; l < gl; l++) */
    /* { */
    /*     plan.gw[l] = g[l - glh]; */
    /* } */

    *pout = plan;
    return status;
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan)** plan)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(*plan);
    LTFAT_NAME(dgtreal_fb_plan)* pp = *plan;
    LTFAT_SAFEFREEALL(pp->sbuf, pp->cbuf, pp->gw, pp->fw);
    LTFAT_FFTW(destroy_plan)(pp->p_small);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_execute)(LTFAT_NAME(dgtreal_fb_plan)* plan,
                               const LTFAT_REAL* f,
                               const ltfatInt L, const ltfatInt W,
                               LTFAT_COMPLEX* cout)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(f); CHECKNULL(cout);
    CHECK(LTFATERR_BADARG, L >= plan->gl && !(L % plan->a) ,
          "L (passed %d) must be positive and divisible by a (passed %d).", L, plan->a);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");

    /*  --------- initial declarations -------------- */
    const ltfatInt a = plan->a;
    const ltfatInt M = plan->M;
    const ltfatInt N = L / a;

    const ltfatInt gl = plan->gl;
    LTFAT_REAL*    sbuf = plan->sbuf;
    LTFAT_COMPLEX* cbuf = plan->cbuf;
    LTFAT_REAL* fw = plan->fw;

    /* These are floor operations. */
    const ltfatInt glh = plan->gl / 2;
    const ltfatInt M2 = M / 2 + 1;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));

    const LTFAT_REAL* fbd;

    /*  ---------- main body ----------- */

    /*----- Handle the first boundary using periodic boundary conditions.*/
    for (ltfatInt n = 0; n < glh_d_a; n++)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = f + L - (glh - n * a) + L * w;
            for (ltfatInt l = 0; l < glh - n * a; l++)
                fw[l]  = fbd[l] * plan->gw[l];

            fbd = f - (glh - n * a) + L * w;
            for (ltfatInt l = glh - n * a; l < gl; l++)
                fw[l]  = fbd[l] * plan->gw[l];

            THE_SUM_REAL
        }
    }

    /* ----- Handle the middle case. --------------------- */
    for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = f + (n * a - glh + L * w);
            for (ltfatInt l = 0; l < gl; l++)
                fw[l]  = fbd[l] * plan->gw[l];

            THE_SUM_REAL
        }

    }

    /* Handle the last boundary using periodic boundary conditions. */
    for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = f + (n * a - glh + L * w);
            for (ltfatInt l = 0; l < L - n * a + glh; l++)
                fw[l]  = fbd[l] * plan->gw[l];

            fbd = f - (L - n * a + glh) + L * w;
            for (ltfatInt l = L - n * a + glh; l < gl; l++)
                fw[l]  = fbd[l] * plan->gw[l];

            THE_SUM_REAL
        }
    }

error:
    return status;
}

#undef THE_SUM_REAL
