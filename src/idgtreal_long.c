#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME(idgtreal_long_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt W;
    ltfatInt c;
    ltfatInt h_a;
    ltfat_phaseconvention ptype;
    LTFAT_REAL scalconst;
    LTFAT_REAL* f;
    LTFAT_COMPLEX* cin;
    LTFAT_COMPLEX* gf;
    LTFAT_COMPLEX* ff;
    LTFAT_COMPLEX* cf;
    LTFAT_COMPLEX* cbuf;
    LTFAT_REAL* cwork;
    LTFAT_REAL* sbuf;
    LTFAT_NAME(ifftreal_plan)* p_veryend;
    LTFAT_NAME(ifftreal_plan)* p_before;
    LTFAT_NAME(fftreal_plan)* p_after;
};


LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long)(const LTFAT_COMPLEX* cin, const LTFAT_REAL* g,
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M,
                          const ltfat_phaseconvention ptype, LTFAT_REAL* f)
{
    LTFAT_NAME(idgtreal_long_plan)* plan = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cin); CHECKNULL(f);

    CHECKSTATUS(
        LTFAT_NAME(idgtreal_long_init)((LTFAT_COMPLEX*)cin, g, L, W, a, M, f,
                                       ptype, FFTW_ESTIMATE, &plan),
        "Init failed");

    LTFAT_NAME(idgtreal_long_execute)(plan);

    LTFAT_NAME(idgtreal_long_done)(&plan);

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_init)(LTFAT_COMPLEX* cin, const LTFAT_REAL* g,
                               const ltfatInt L, const ltfatInt W,
                               const ltfatInt a, const ltfatInt M, LTFAT_REAL* f,
                               const ltfat_phaseconvention ptype, unsigned flags,
                               LTFAT_NAME(idgtreal_long_plan)** pout)
{
    ltfatInt minL, h_m, b, N, p, q, d, d2, size;
    // Downcast to int
    LTFAT_NAME(idgtreal_long_plan)* plan = NULL;
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NULLPOINTER, (flags & FFTW_ESTIMATE) || cin != NULL,
          "cin cannot be NULL if flags is not FFTW_ESTIMATE");

    // CHECKNULL(f); // can be NULL
    CHECKNULL(g); CHECKNULL(pout);
    CHECK(LTFATERR_BADSIZE, L > 0, "L (passed %d) must be positive", L);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W (passed %d) must be positive.", W);
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a (passed %d) must be positive.", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M (passed %d) must be positive.", M);
    CHECK(LTFATERR_CANNOTHAPPEN, ltfat_phaseconvention_is_valid(ptype),
          "Invalid ltfat_phaseconvention enum value." );

    minL = ltfat_lcm(a, M);
    CHECK(LTFATERR_BADTRALEN,
          !(L % minL), "L must be divisible by lcm(a,M)=%d.", minL);

    CHECKMEM(plan =
                 (LTFAT_NAME(idgtreal_long_plan)*) ltfat_calloc(1, sizeof * plan));

    /*  ----------- calculation of parameters and plans -------- */

    plan->a = a; plan->L = L; plan->M = M; plan->W = W; plan->ptype = ptype;
    b = L / M;
    N = L / a;

    plan->c = ltfat_gcd(a, M, &plan->h_a, &h_m);
    p = a / plan->c;
    q = M / plan->c;
    d = b / p;

    /* This is a floor operation. */
    d2 = d / 2 + 1;

    size = wfacreal_size(L, a, M);
    CHECKMEM( plan->gf    = LTFAT_NAME_COMPLEX(malloc)(size));
    CHECKMEM( plan->ff    = LTFAT_NAME_COMPLEX(malloc)(d2 * p * q * W));
    CHECKMEM( plan->cf    = LTFAT_NAME_COMPLEX(malloc)(d2 * q * q * W));
    CHECKMEM( plan->cwork = LTFAT_NAME_REAL(malloc)(M * N * W));
    CHECKMEM( plan->cbuf  = LTFAT_NAME_COMPLEX(malloc)(d2));
    CHECKMEM( plan->sbuf  = LTFAT_NAME_REAL(malloc)(d));
    plan->cin = cin;
    plan->f = f;

    LTFAT_NAME(wfacreal)(g, L, 1, a, M, plan->gf);

    /* Scaling constant needed because of FFTWs normalization. */
    plan->scalconst = (LTFAT_REAL) ( 1.0 / ((double)d * sqrt((double)M)));

    CHECKSTATUS(
        LTFAT_NAME(ifftreal_init)(d, 1, plan->cbuf, plan->sbuf, flags, &plan->p_before),
        "FFTW plan failed.");

    CHECKSTATUS(
        LTFAT_NAME(fftreal_init)(d, 1, plan->sbuf, plan->cbuf, flags, &plan->p_after),
        "FFTW plan failed.");

    CHECKSTATUS(
        LTFAT_NAME(ifftreal_init)(M, N * W, plan->cin, plan->cwork,
                                  flags | FFTW_PRESERVE_INPUT, &plan->p_veryend),
        "FFTW plan failed.");

    *pout = plan;
    return status;
error:
    if (plan)
    {
        if (plan->p_before) LTFAT_NAME(ifftreal_done)(&plan->p_before);
        if (plan->p_after) LTFAT_NAME(fftreal_done)(&plan->p_after);
        if (plan->p_veryend) LTFAT_NAME(ifftreal_done)(&plan->p_veryend);
        LTFAT_SAFEFREEALL(plan->gf, plan->ff, plan->cf, plan->cwork, plan->cbuf,
                          plan->sbuf);
        ltfat_free(plan);
    }
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_execute)(LTFAT_NAME(idgtreal_long_plan)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(p->f); CHECKNULL(p->cin);

    LTFAT_NAME(ifftreal_execute)(p->p_veryend);

    if (p->ptype)
        LTFAT_NAME_REAL(dgtphaseunlockhelper)(p->cwork, p->L, p->W, p->a, p->M,
                                              p->cwork);

    LTFAT_NAME(idgtreal_walnut_execute)(p);
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_execute_newarray)(LTFAT_NAME(idgtreal_long_plan)* p,
        const LTFAT_COMPLEX* c, LTFAT_REAL* f)
{
    LTFAT_NAME(idgtreal_long_plan) p2;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(c); CHECKNULL(f);

    // The plan was created with the FFTW_PRESERVE_INPUT so it is ok to cast away the const
    LTFAT_NAME(ifftreal_execute_newarray)(p->p_veryend, c, p->cwork);

    if (p->ptype)
        LTFAT_NAME_REAL(dgtphaseunlockhelper)(p->cwork, p->L, p->W, p->a, p->M,
                                              p->cwork);

    // Make a shallow copy and rewrite f
    p2 = *p;
    p2.f = f;

    LTFAT_NAME(idgtreal_walnut_execute)(&p2);
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_done)(LTFAT_NAME(idgtreal_long_plan)** plan)
{
    LTFAT_NAME(idgtreal_long_plan)* p;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(*plan);
    p = *plan;

    LTFAT_NAME(ifftreal_done)(&p->p_before);
    LTFAT_NAME(fftreal_done)(&p->p_after);
    LTFAT_NAME(ifftreal_done)(&p->p_veryend);
    LTFAT_SAFEFREEALL(p->gf, p->ff, p->cf, p->cwork, p->cbuf, p->sbuf);

    ltfat_free(p);
    p = NULL;
error:
    return status;
}

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_walnut_execute)(LTFAT_NAME(idgtreal_long_plan)* p)
{
    const ltfatInt b = p->L / p->M;
    const ltfatInt N = p->L / p->a;

    const ltfatInt c = p->c;
    const ltfatInt pp = p->a / c;
    const ltfatInt q = p->M / c;
    const ltfatInt d = b / pp;

    const ltfatInt L = p->L;
    const ltfatInt W = p->W;
    const ltfatInt a = p->a;
    const ltfatInt M = p->M;

    /* This is a floor operation. */
    const ltfatInt d2 = d / 2 + 1;

    ltfatInt h_a = -p->h_a;

    /* Scaling constant needed because of FFTWs normalization. */
    const LTFAT_REAL scalconst = p->scalconst;

    const ltfatInt ld4c = p->M * N;

    /* Leading dimensions of cf */
    const ltfatInt ld3b = q * q * W;

    /* Leading dimensions of the 4dim array. */
    const ltfatInt ld2ff = pp * q * W;

    /* -------- Main loop ----------------------------------- */
    for (ltfatInt r = 0; r < c; r++)
    {
        /* -------- compute coefficient factorization ----------- */

        LTFAT_COMPLEX* cfp = p->cf;

        for (ltfatInt w = 0; w < W; w++)
        {
            /* Complete inverse fac of coefficients */
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt u = 0; u < q; u++)
                {
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        p->sbuf[s] = p->cwork[r + l * c +
                                              ltfat_positiverem(u + s * q - l * h_a, N) *
                                              M + w * ld4c];
                    }

                    /* Do inverse fft of length d */
                    LTFAT_NAME(fftreal_execute)(p->p_after);

                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        cfp[s * ld3b] = p->cbuf[s];
                    }
                    /* Advance the cf pointer. This is only done in this
                    * one place, because the loops are placed such that
                    * this pointer will advance linearly through
                    * memory. Reordering the loops will break this. */
                    cfp++;
                }
            }
        }
        /* -------- compute matrix multiplication ---------- */
        /* Do the matmul  */
        for (ltfatInt s = 0; s < d2; s++)
        {
            const LTFAT_COMPLEX* gbase = p->gf + (r + s * c) * pp * q;
            LTFAT_COMPLEX*       fbase = p->ff + s * pp * q * W;
            const LTFAT_COMPLEX* cbase = p->cf + s * q * q * W;

            for (ltfatInt nm = 0; nm < q * W; nm++)
            {
                for (ltfatInt km = 0; km < pp; km++)
                {
                    fbase[km + nm * pp] = 0.0;
                    for (ltfatInt mm = 0; mm < q; mm++)
                    {
                        fbase[km + nm * pp] += gbase[km + mm * pp] * cbase[mm + nm * q];
                    }
                    /* Scale because of FFTWs normalization. */
                    fbase[km + nm * pp] = fbase[km + nm * pp] * scalconst;
                }
            }
        }
        /* ----------- compute inverse signal factorization ---------- */

        LTFAT_COMPLEX* ffp = p->ff;
        LTFAT_REAL*    fp  = p->f + r;

        for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt k = 0; k < pp; k++)
                {
                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        p->cbuf[s] = ffp[s * ld2ff];
                    }

                    LTFAT_NAME(ifftreal_execute)(p->p_before);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        fp[ltfat_positiverem(k * M + s * pp * M - l * h_a * a, L)] = p->sbuf[s];
                    }

                    /* Advance the ff pointer. This is only done in this
                    * one place, because the loops are placed such that
                    * this pointer will advance linearly through
                    * memory. Reordering the loops will break this. */
                    ffp++;
                }
            }
            fp += L;
        }
    }
}
