#include "dgt_long_private.h"

LTFAT_EXTERN int
LTFAT_NAME(dgt_long)(const LTFAT_TYPE* f, const LTFAT_TYPE* g,
                     const ltfatInt L, const ltfatInt W,
                     const ltfatInt a, const ltfatInt M,
                     const ltfat_phaseconvention ptype, LTFAT_COMPLEX* cout)
{
    LTFAT_NAME(dgt_long_plan)* plan = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        LTFAT_NAME(dgt_long_init)(f, g, L, W, a, M, cout, ptype, FFTW_ESTIMATE, &plan)
        , "Init failed");

    // Nothing should go wrong once the plan is created
    LTFAT_NAME(dgt_long_execute)(plan);

    LTFAT_NAME(dgt_long_done)(&plan);
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(dgt_long_init)(const LTFAT_TYPE* f, const LTFAT_TYPE* g,
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX* cout,
                          const ltfat_phaseconvention ptype, unsigned flags,
                          LTFAT_NAME(dgt_long_plan)** pout)
{
    LTFAT_NAME(dgt_long_plan)* plan = NULL;
    ltfatInt h_m, N, b, p, q, d, minL;

    int status = LTFATERR_SUCCESS;
    // CHECKNULL(f); // Can be NULL
    CHECK(LTFATERR_NULLPOINTER, (flags & FFTW_ESTIMATE) || cout != NULL,
          "cout cannot be NULL if flags is not FFTW_ESTIMATE");

    CHECKNULL(g); CHECKNULL(pout);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W (passed %d) must be positive.", W);
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a (passed %d) must be positive.", a);
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M (passed %d) must be positive.", M);

    minL = ltfat_lcm(a, M);
    CHECK(LTFATERR_BADARG,
          L > 0  && !(L % minL),
          "L (passed %d) must be positive and divisible by lcm(a,M)=%d.",
          L, minL);

    CHECKMEM(plan = (LTFAT_NAME(dgt_long_plan)*)ltfat_calloc(1, sizeof * plan));

    plan->a = a;
    plan->M = M;
    plan->L = L;
    plan->W = W;
    plan->ptype = ptype;
    N = L / a;
    b = L / M;

    plan->c = ltfat_gcd(a, M, &plan->h_a, &h_m);
    p = a / plan->c;
    q = M / plan->c;
    d = b / p;
    plan->h_a = -plan->h_a;

    CHECKMEM( plan->sbuf = LTFAT_NAME_REAL(malloc)(2 * d));
    CHECKMEM( plan->gf   = LTFAT_NAME_COMPLEX(malloc)(L));
    CHECKMEM( plan->ff = LTFAT_NAME_REAL(malloc)(2 * d * p * q * W));
    CHECKMEM( plan->cf = LTFAT_NAME_REAL(malloc)(2 * d * q * q * W));
    plan->cout = cout;
    plan->f    = f;

    /* Get factorization of window */
    CHECKSTATUS(
        LTFAT_NAME(wfac)(g, L, 1, a, M, plan->gf),
        "wfac call failed.");

    CHECKSTATUS(
        LTFAT_NAME_REAL(fft_init)(M, N * W, cout, cout, flags, &plan->p_veryend),
        "FFTW plan creation failed.");

    CHECKSTATUS(
        LTFAT_NAME_REAL(fft_init)(d, 1, (LTFAT_COMPLEX*) plan->sbuf,
                                  (LTFAT_COMPLEX*) plan->sbuf, flags, &plan->p_before),
        "FFTW plan creation failed.");

    CHECKSTATUS(
        LTFAT_NAME_REAL(ifft_init)(d, 1, (LTFAT_COMPLEX*) plan->sbuf,
                                   (LTFAT_COMPLEX*) plan->sbuf, flags, &plan->p_after),
        "FFTW plan creation failed.");

    // Assign the "return" value
    *pout = plan;
    return status;
error:
    if (plan)
    {
        if (plan->p_veryend) LTFAT_NAME_REAL(fft_done)(&plan->p_veryend);
        if (plan->p_before)  LTFAT_NAME_REAL(fft_done)(&plan->p_before);
        if (plan->p_after)   LTFAT_NAME_REAL(ifft_done)(&plan->p_after);
        LTFAT_SAFEFREEALL(plan->sbuf, plan->gf, plan->ff, plan->cf);
        ltfat_free(plan);
    }
    *pout = NULL;
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(dgt_long_execute)(LTFAT_NAME(dgt_long_plan)* plan)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(plan->f); CHECKNULL(plan->cout);

    LTFAT_NAME(dgt_walnut_execute)(plan, plan->cout);

    if (LTFAT_TIMEINV == plan->ptype)
        LTFAT_NAME_COMPLEX(dgtphaselockhelper)(plan->cout, plan->L, plan->W,
                                               plan->a, plan->M, plan->cout);

    /* FFT to modulate the coefficients. */
    LTFAT_NAME_REAL(fft_execute)(plan->p_veryend);

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(dgt_long_execute_newarray)(LTFAT_NAME(dgt_long_plan)* plan,
                                      const LTFAT_TYPE f[], LTFAT_COMPLEX c[])
{
    LTFAT_NAME(dgt_long_plan) plan2;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(f); CHECKNULL(c);

    // Make a shallow copy and assign f
    plan2 = *plan;
    plan2.f = f;

    LTFAT_NAME(dgt_walnut_execute)(&plan2, c);

    if (LTFAT_TIMEINV == plan->ptype)
        LTFAT_NAME_COMPLEX(dgtphaselockhelper)(c, plan->L, plan->W,
                                               plan->a, plan->M, c);

    /* FFT to modulate the coefficients. */
    LTFAT_NAME_REAL(fft_execute_newarray)(plan->p_veryend, c, c);

    /* LTFAT_FFTW(execute_dft)(plan->p_veryend, (LTFAT_FFTW(complex)*)c, */
    /*                         (LTFAT_FFTW(complex)*)c); */

error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan)** plan)
{
    LTFAT_NAME(dgt_long_plan)* pp;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(*plan);
    pp = *plan;

    LTFAT_NAME_REAL(fft_done)(&pp->p_veryend);
    LTFAT_NAME_REAL(fft_done)(&pp->p_before);
    LTFAT_NAME_REAL(ifft_done)(&pp->p_after);
    LTFAT_SAFEFREEALL(pp->sbuf, pp->gf, pp->ff, pp->cf);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}



/*  This routine computes the DGT factorization using strided FFTs so
    the memory layout is optimized for the matrix product. Compared to
    dgt_fac_1, it moves the r-loop to be the outermost loop to
    conserve memory and hopefully use the cache hierachy better

    The routine uses a very small buffer to do the DFTs.

    Integer indexing is optimized.

    Special code for integer oversampling.

    Code works on LTFAT_REAL's instead on LTFAT_COMPLEX
*/


LTFAT_EXTERN int
LTFAT_NAME(dgt_walnut_execute)(LTFAT_NAME(dgt_long_plan)* plan,
                               LTFAT_COMPLEX* cout)
{

    /*  --------- initial declarations -------------- */

    LTFAT_REAL* gbase, *fbase, *cbase;

    ltfatInt rem;

    LTFAT_REAL* ffp, *cfp;
    LTFAT_TYPE* fp;

    /*  ----------- calculation of parameters and plans -------- */

    const ltfatInt a = plan->a;
    const ltfatInt M = plan->M;
    const ltfatInt L = plan->L;
    const ltfatInt W = plan->W;
    const ltfatInt N = L / a;
    const ltfatInt c = plan->c;
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = N / q;

    LTFAT_TYPE* f = (LTFAT_TYPE*) plan->f;
    const LTFAT_COMPLEX* gf = (const LTFAT_COMPLEX*)plan->gf;

    const ltfatInt h_a = plan->h_a;

    LTFAT_REAL* sbuf = plan->sbuf;
    //LTFAT_COMPLEX* cout = plan->cout;

    /* Scaling constant needed because of FFTWs normalization. */
    const LTFAT_REAL scalconst = (const LTFAT_REAL)( 1.0 / ((double)d * sqrt((
                                     double)M)) );

    /* Leading dimensions of the 4dim array. */
    const ltfatInt ld2a = 2 * p * q * W;

    /* Leading dimensions of cf */
    const ltfatInt ld3b = 2 * q * q * W;
    const ltfatInt ld5c = M * N;

    /* --------- main loop begins here ------------------- */
    for (ltfatInt r = 0; r < c; r++)
    {
        /*  ---------- compute signal factorization ----------- */
        ffp = plan->ff;
        fp = f + r;
        if (p == 1)
        {
            /* Integer oversampling case */
            for (ltfatInt w = 0; w < W; w++)
            {
                for (ltfatInt l = 0; l < q; l++)
                {
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        rem = (s * M + l * a) % L;
#ifdef LTFAT_COMPLEXTYPE
                        sbuf[2 * s]   = ltfat_real(fp[rem]);
                        sbuf[2 * s + 1] = ltfat_imag(fp[rem]);
#else
                        sbuf[2 * s]   = fp[rem];
                        sbuf[2 * s + 1] = 0.0;
#endif
                    }

                    LTFAT_NAME_REAL(fft_execute)(plan->p_before);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        ffp[s * ld2a]   = sbuf[2 * s] * scalconst;
                        ffp[s * ld2a + 1] = sbuf[2 * s + 1] * scalconst;
                    }
                    ffp += 2;
                }
                fp += L;
            }
            fp -= L * W;

            /* Do the Matmul */
            for (ltfatInt s = 0; s < d; s++)
            {
                gbase = (LTFAT_REAL*)gf + 2 * (r + s * c) * q;
                fbase = plan->ff + 2 * s * q * W;
                cbase = plan->cf + 2 * s * q * q * W;

                for (ltfatInt nm = 0; nm < q * W; nm++)
                {
                    for (ltfatInt mm = 0; mm < q; mm++)
                    {
                        cbase[0] = gbase[0] * fbase[0] + gbase[1] * fbase[1];
                        cbase[1] = gbase[0] * fbase[1] - gbase[1] * fbase[0];
                        gbase += 2;
                        cbase += 2;
                    }
                    gbase -= 2 * q;
                    fbase += 2;
                }
                cbase -= 2 * q * q * W;
            }
        }
        else
        {
            /* rational sampling case */

            for (ltfatInt w = 0; w < W; w++)
            {
                for (ltfatInt l = 0; l < q; l++)
                {
                    for (ltfatInt k = 0; k < p; k++)
                    {
                        for (ltfatInt s = 0; s < d; s++)
                        {
                            rem = ltfat_positiverem(k * M + s * p * M - l * h_a * a, L);
#ifdef LTFAT_COMPLEXTYPE
                            sbuf[2 * s]   = ltfat_real(fp[rem]);
                            sbuf[2 * s + 1] = ltfat_imag(fp[rem]);
#else
                            sbuf[2 * s]   = fp[rem];
                            sbuf[2 * s + 1] = 0.0;
#endif
                        }

                        LTFAT_NAME_REAL(fft_execute)(plan->p_before);

                        for (ltfatInt s = 0; s < d; s++)
                        {
                            ffp[s * ld2a]   = sbuf[2 * s] * scalconst;
                            ffp[s * ld2a + 1] = sbuf[2 * s + 1] * scalconst;
                        }
                        ffp += 2;
                    }
                }
                fp += L;
            }
            fp -= L * W;

            // Matmul
            for (ltfatInt s = 0; s < d; s++)
            {
                gbase = (LTFAT_REAL*)gf + 2 * (r + s * c) * p * q;
                fbase = plan->ff + 2 * s * p * q * W;
                cbase = plan->cf + 2 * s * q * q * W;

                for (ltfatInt nm = 0; nm < q * W; nm++)
                {
                    for (ltfatInt mm = 0; mm < q; mm++)
                    {
                        cbase[0] = 0.0;
                        cbase[1] = 0.0;
                        for (ltfatInt km = 0; km < p; km++)
                        {
                            cbase[0] += gbase[0] * fbase[0] + gbase[1] * fbase[1];
                            cbase[1] += gbase[0] * fbase[1] - gbase[1] * fbase[0];
                            gbase += 2;
                            fbase += 2;
                        }
                        fbase -= 2 * p;
                        cbase += 2;
                    }
                    gbase -= 2 * q * p;
                    fbase += 2 * p;
                }
                cbase -= 2 * q * q * W;
                fbase -= 2 * p * q * W;
            }

        } /* end of if p==1 */

        /*  -------  compute inverse coefficient factorization ------- */
        cfp = plan->cf;

        /* Cover both integer and rational sampling case */
        for (ltfatInt w = 0; w < W; w++)
        {
            /* Complete inverse fac of coefficients */
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt u = 0; u < q; u++)
                {
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        sbuf[2 * s]   = cfp[s * ld3b];
                        sbuf[2 * s + 1] = cfp[s * ld3b + 1];
                    }
                    cfp += 2;

                    /* Do inverse fft of length d */
                    LTFAT_NAME_REAL(ifft_execute)(plan->p_after);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        rem = r + l * c + ltfat_positiverem(u + s * q - l * h_a, N) * M + w * ld5c;
                        LTFAT_REAL* coutTmp = (LTFAT_REAL*) &cout[rem];
                        coutTmp[0] = sbuf[2 * s];
                        coutTmp[1] = sbuf[2 * s + 1];
                    }
                }
            }
        }


        /* ----------- Main loop ends here ------------------------ */
    }

    return LTFATERR_SUCCESS;
}
