#include "dgtreal_long_private.h"

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long)(const LTFAT_REAL* f, const LTFAT_REAL* g,
                         const ltfatInt L, const ltfatInt W, const ltfatInt a,
                         const ltfatInt M, const dgt_phasetype ptype,
                         LTFAT_COMPLEX* cout)
{

    LTFAT_NAME(dgtreal_long_plan) *plan = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKSTATUS(
        LTFAT_NAME(dgtreal_long_init)(f, g, L, W, a, M, cout, ptype, FFTW_ESTIMATE,
                                      &plan),
        "Init failed");

    LTFAT_NAME(dgtreal_long_execute)(plan);

    LTFAT_NAME(dgtreal_long_done)(&plan);

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_init)(const LTFAT_REAL* f, const LTFAT_REAL* g,
                              const ltfatInt L, const ltfatInt W,
                              const ltfatInt a, const ltfatInt M,
                              LTFAT_COMPLEX* cout, const dgt_phasetype ptype,
                              unsigned flags, LTFAT_NAME(dgtreal_long_plan)** pout)
{
    LTFAT_NAME(dgtreal_long_plan)* plan = NULL;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(f); CHECKNULL(g); CHECKNULL(cout); CHECKNULL(pout);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    ltfatInt minL = ltfat_lcm(a, M);
    CHECK(LTFATERR_BADARG,
          L > 0  && !(L % minL),
          "L (passed %d) must be positive and divisible by lcm(a,M)=%d.",
          L, minL);

    CHECKMEM(plan = ltfat_calloc(1, sizeof * plan));

    ltfatInt h_m;

    plan->a = a;
    plan->M = M;
    plan->L = L;
    plan->W = W;
    plan->ptype = ptype;
    const ltfatInt N = L / a;

    plan->c = ltfat_gcd(a, M, &plan->h_a, &h_m);
    const ltfatInt b = L / M;
    const ltfatInt p = a / plan->c;
    const ltfatInt q = M / plan->c;
    const ltfatInt d = b / p;
    plan->h_a = -plan->h_a;

    const ltfatInt M2 = M / 2 + 1;
    const ltfatInt d2 = d / 2 + 1;
    const ltfatInt wfs = wfacreal_size(L, a, M);

    plan->cout = cout;
    plan->f    = f;
    CHECKMEM( plan->sbuf = ltfat_malloc( d * sizeof * plan->sbuf ));
    CHECKMEM( plan->cbuf = ltfat_malloc(d2 * sizeof * plan->cbuf ));
    CHECKMEM( plan->ff = ltfat_malloc(2 * d2 * p * q * W * sizeof * plan->ff));
    CHECKMEM( plan->cf = ltfat_malloc(2 * d2 * q * q * W * sizeof * plan->cf));
    CHECKMEM( plan->gf = ltfat_malloc(wfs * sizeof * plan->gf));
    CHECKMEM( plan->cwork = ltfat_malloc(M * N * W * sizeof * plan->cwork));

    /* Get factorization of window */
    LTFAT_NAME(wfacreal)(g, L, 1, a, M, plan->gf);

    /* Create plans. In-place. */
    // Downcast to int
    int Mint = (int) plan->M;

    plan->p_veryend =
        LTFAT_FFTW(plan_many_dft_r2c)(1, &Mint, N * W,
                                      plan->cwork, NULL, 1, M,
                                      cout, NULL, 1, M2, flags);

    CHECKINIT(plan->p_veryend, "FFTW plan creation failed." );

    plan->p_before =
        LTFAT_FFTW(plan_dft_r2c_1d)(d, plan->sbuf, plan->cbuf, flags);

    CHECKINIT(plan->p_before, "FFTW plan creation failed." );

    plan->p_after  =
        LTFAT_FFTW(plan_dft_c2r_1d)(d, plan->cbuf, plan->sbuf, flags);

    CHECKINIT(plan->p_after, "FFTW plan creation failed.");

    *pout = plan;
    return status;
error:
    if (plan)
    {
        if (plan->p_veryend) LTFAT_FFTW(destroy_plan)(plan->p_veryend);
        if (plan->p_before)  LTFAT_FFTW(destroy_plan)(plan->p_before);
        if (plan->p_after)   LTFAT_FFTW(destroy_plan)(plan->p_after);
        LTFAT_SAFEFREEALL(plan->sbuf, plan->cbuf, plan->gf, plan->ff, plan->cf,
                          plan->cwork);
        ltfat_free(plan);
    }
    return status;
}




LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_execute)(LTFAT_NAME(dgtreal_long_plan)* plan)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan);

    LTFAT_NAME(dgtreal_walnut_plan)(plan);

    if (plan->ptype)
        LTFAT_NAME_REAL(dgtphaselockhelper)(plan->cwork, plan->L, plan->W, plan->a,
                                            plan->M, plan->cwork);

    /* FFT to modulate the coefficients. */
    LTFAT_FFTW(execute)(plan->p_veryend);

error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan)** plan)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(plan); CHECKNULL(*plan);
    LTFAT_NAME(dgtreal_long_plan)* pp = *plan;
    LTFAT_FFTW(destroy_plan)(pp->p_veryend);
    LTFAT_FFTW(destroy_plan)(pp->p_before);
    LTFAT_FFTW(destroy_plan)(pp->p_after);
    LTFAT_SAFEFREEALL(pp->sbuf, pp->cbuf, pp->cwork,
                      pp->gf, pp->ff, pp->cf);
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

*/

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan)* plan)
{
    /*  --------- initial declarations -------------- */

    const ltfatInt a = plan->a;
    const ltfatInt M = plan->M;
    const ltfatInt L = plan->L;
    const ltfatInt W = plan->W;
    const ltfatInt N = L / a;
    const ltfatInt c = plan->c;
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = N / q;


    /* This is a floor operation. */
    const ltfatInt d2 = d / 2 + 1;

    const LTFAT_REAL* f = plan->f;
    const LTFAT_COMPLEX* gf = (const LTFAT_COMPLEX*)plan->gf;

    const ltfatInt h_a = plan->h_a;

    LTFAT_REAL* sbuf = plan->sbuf;
    LTFAT_COMPLEX* cbuf = plan->cbuf;

    LTFAT_REAL* cout = plan->cwork;


    LTFAT_REAL* gbase, *fbase, *cbase;

    LTFAT_REAL* ffp;

    const LTFAT_REAL* fp;

    /* Scaling constant needed because of FFTWs normalization. */
    const LTFAT_REAL scalconst = 1.0 / ((LTFAT_REAL)d * sqrt((LTFAT_REAL)M));

    /* Leading dimensions of the 4dim array. */
    const ltfatInt ld2a = 2 * p * q * W;

    /* Leading dimensions of cf */
    const ltfatInt ld3b = 2 * q * q * W;

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
                        sbuf[s]   = fp[(s * M + l * a) % L];
                    }

                    LTFAT_FFTW(execute)(plan->p_before);

                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        ffp[s * ld2a]   = LTFAT_COMPLEXH(creal)(cbuf[s]) * scalconst;
                        ffp[s * ld2a + 1] = LTFAT_COMPLEXH(cimag)(cbuf[s]) * scalconst;
                    }
                    ffp += 2;
                }
                fp += L;
            }
            /* fp -= 2 * L * W; */
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
                            sbuf[s]   = fp[ positiverem(k * M + s * p * M - l * h_a * a, L) ];
                        }

                        LTFAT_FFTW(execute)(plan->p_before);

                        for (ltfatInt s = 0; s < d2; s++)
                        {
                            ffp[s * ld2a]   = LTFAT_COMPLEXH(creal)(cbuf[s]) * scalconst;
                            ffp[s * ld2a + 1] = LTFAT_COMPLEXH(cimag)(cbuf[s]) * scalconst;
                        }
                        ffp += 2;
                    }
                }
                fp += L;
            }
            /* fp -= 2 * L * W; */
        }

        /* ----------- compute matrix multiplication ----------- */

        /* Do the matmul  */
        if (p == 1)
        {
            /* Integer oversampling case */


            /* Rational oversampling case */
            for (ltfatInt s = 0; s < d2; s++)
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

            /* Rational oversampling case */
            for (ltfatInt s = 0; s < d2; s++)
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
        }



        /*  -------  compute inverse coefficient factorization ------- */
        LTFAT_REAL* cfp = plan->cf;
        const ltfatInt ld5c = M * N;

        /* Cover both integer and rational sampling case */
        for (ltfatInt w = 0; w < W; w++)
        {
            /* Complete inverse fac of coefficients */
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt u = 0; u < q; u++)
                {
                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        LTFAT_REAL* cbufTmp = (LTFAT_REAL*) &cbuf[s];
                        cbufTmp[0] = cfp[s * ld3b];
                        cbufTmp[1] = cfp[s * ld3b + 1];
                    }
                    cfp += 2;

                    /* Do inverse fft of length d */
                    LTFAT_FFTW(execute)(plan->p_after);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        cout[ r + l * c + positiverem(u + s * q - l * h_a, N)*M + w * ld5c ] = sbuf[s];
                    }
                }
            }
        }


        /* ----------- Main loop ends here ------------------------ */
    }

    return 0;
}
