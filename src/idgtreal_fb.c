#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME(idgtreal_fb_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt gl;
    ltfat_phaseconvention ptype;
    LTFAT_COMPLEX* cbuf;
    LTFAT_REAL*    crbuf;
    LTFAT_REAL*    gw;
    LTFAT_REAL*    ff;
    LTFAT_FFTW(plan) p_small;
};


/* ------------------- IDGTREAL ---------------------- */

#define THE_SUM_REAL { \
    memcpy(cbuf,cin+n*M2+w*M2*N, M2*sizeof*cbuf); \
    LTFAT_FFTW(execute)(p->p_small); \
    LTFAT_NAME_REAL(circshift)(crbuf,M,p->ptype?glh:-n*a+glh,ff); \
    LTFAT_NAME_REAL(periodize_array)(ff,M,gl,ff); \
    for (ltfatInt ii=0; ii<gl; ii++) \
        ff[ii] *= gw[ii]; \
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX* cin, const LTFAT_REAL* g,
                        const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M,
                        const ltfat_phaseconvention ptype, LTFAT_REAL* f)

{
    LTFAT_NAME(idgtreal_fb_plan)* plan = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        LTFAT_NAME(idgtreal_fb_init)(g, gl, a, M, ptype, FFTW_ESTIMATE, &plan),
        "Init failed");

    CHECKSTATUS(
        LTFAT_NAME(idgtreal_fb_execute)(plan, cin, L, W, f),
        "Execute failed");

error:
    if (plan) LTFAT_NAME(idgtreal_fb_done)(&plan);
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_init)(const LTFAT_REAL* g, const ltfatInt gl,
                             const ltfatInt a, const ltfatInt M, const ltfat_phaseconvention ptype,
                             unsigned flags, LTFAT_NAME(idgtreal_fb_plan)** pout)
{
    LTFAT_NAME(idgtreal_fb_plan)* p = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g); CHECKNULL(pout);
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl (passed %d) must be positive.", gl);
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a (passed %d) must be positive.", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M (passed %d) must be positive.", M);

    CHECKMEM(p = ltfat_calloc(1, sizeof * p));

    p->ptype = ptype;
    p->a = a;
    p->M = M;
    p->gl = gl;

    /* This is a floor operation. */
    const ltfatInt M2 = M / 2 + 1;

    CHECKMEM( p->cbuf  = ltfat_malloc(M2 * sizeof * p->cbuf));
    CHECKMEM( p->crbuf = ltfat_malloc( M * sizeof * p->crbuf));
    CHECKMEM( p->gw    = ltfat_malloc(gl * sizeof * p->gw));
    CHECKMEM( p->ff    = ltfat_malloc((gl > M ? gl : M) * sizeof * p->ff));

    /* Create plan. In-place. */
    p->p_small = LTFAT_FFTW(plan_dft_c2r_1d)(M, p->cbuf, p->crbuf,
                 flags);

    CHECKINIT(p->p_small, "FFTW plan failed.");

    LTFAT_NAME_REAL(fftshift)(g, gl, p->gw);

    *pout = p;
    return status;
error:
    if (p)
    {
        LTFAT_SAFEFREEALL(p->cbuf, p->crbuf, p->gw, p->ff);
        LTFAT_FFTW(destroy_plan)(p->p_small);
        ltfat_free(p);
    }
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_execute)(LTFAT_NAME(idgtreal_fb_plan)* p,
                                const LTFAT_COMPLEX* cin,
                                const ltfatInt L, const ltfatInt W, LTFAT_REAL* f)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(cin); CHECKNULL(f);
    CHECK(LTFATERR_BADARG, L >= p->gl && !(L % p->a) ,
          "L (passed %d) must be positive and divisible by a (passed %d).", L, p->a);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W (passed %d) must be positive.", W);

    const ltfatInt M = p->M;
    const ltfatInt a = p->a;
    const ltfatInt gl = p->gl;
    const ltfatInt N = L / a;

    /* This is a floor operation. */
    const ltfatInt M2 = M / 2 + 1;

    ltfatInt ep, sp;

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));

    LTFAT_COMPLEX* cbuf  = p->cbuf;
    LTFAT_REAL*    crbuf = p->crbuf;
    LTFAT_REAL* gw  = p->gw;
    LTFAT_REAL* ff  = p->ff;

    memset(f, 0, L * W * sizeof * f);

    for (ltfatInt w = 0; w < W; w++)
    {
        LTFAT_REAL* fw = f + w * L;

        /* ----- Handle the first boundary using periodic boundary conditions. --- */
        for (ltfatInt n = 0; n < glh_d_a; n++)
        {
            THE_SUM_REAL;

            sp = ltfat_positiverem(n * a - glh, L);
            ep = ltfat_positiverem(n * a - glh + gl - 1, L);

            /* % Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
                fw[sp + ii] += ff[ii];

            for (ltfatInt ii = 0; ii < ep + 1; ii++)
                fw[ii] += ff[L - sp + ii];
        }


        /* ----- Handle the middle case. --------------------- */
        for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
        {
            THE_SUM_REAL;

            sp = ltfat_positiverem(n * a - glh, L);
            ep = ltfat_positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < ep - sp + 1; ii++)
                fw[ii + sp] += ff[ii];
        }

        /* Handle the last boundary using periodic boundary conditions. */
        for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
        {
            THE_SUM_REAL;

            sp = ltfat_positiverem(n * a - glh, L);
            ep = ltfat_positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
                fw[sp + ii] += ff[ii];

            for (ltfatInt ii = 0; ii < ep + 1; ii++)
                fw[ii] += ff[L - sp + ii];
        }
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_done)(LTFAT_NAME(idgtreal_fb_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    LTFAT_NAME(idgtreal_fb_plan)* pp = *p;
    LTFAT_SAFEFREEALL(pp->cbuf, pp->crbuf, pp->ff, pp->gw);
    LTFAT_FFTW(destroy_plan)(pp->p_small);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}


#undef THE_SUM_REAL
