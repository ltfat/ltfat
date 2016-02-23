#include "ltfat.h"
#include "ltfat_types.h"




LTFAT_EXTERN void
LTFAT_NAME(idgt_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf,
                     const ltfatInt L, const ltfatInt W,
                     const ltfatInt a, const ltfatInt M,
                     const dgt_phasetype ptype,
                     LTFAT_COMPLEX *f)
{

    /*  --------- initial declarations -------------- */

    ltfatInt h_a, h_m;

    LTFAT_FFTW(plan) p_before, p_after, p_veryend;
    LTFAT_COMPLEX *ff, *cf, *cwork, *cbuf;

    /*  ----------- calculation of parameters and plans -------- */

    const ltfatInt b = L / M;
    const ltfatInt N = L / a;

    const ltfatInt c = gcd(a, M, &h_a, &h_m);
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = b / p;

    h_a = -h_a;

    ff    = ltfat_malloc(d * p * q * W * sizeof * ff);
    cf    = ltfat_malloc(d * q * q * W * sizeof * cf);
    cwork = ltfat_malloc(M * N * W * sizeof * cwork);
    cbuf  = ltfat_malloc(d * sizeof * cbuf);

    /* Scaling constant needed because of FFTWs normalization. */
    const LTFAT_REAL scalconst = 1.0 / ((LTFAT_REAL)d * sqrt((LTFAT_REAL)M));

    /* Create plans. In-place. */

    p_after  = LTFAT_FFTW(plan_dft_1d)(d, cbuf, cbuf,
                                       FFTW_FORWARD, FFTW_ESTIMATE);

    p_before = LTFAT_FFTW(plan_dft_1d)(d, cbuf, cbuf,
                                       FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create plan. Copy data so we do not overwrite input. Therefore
       it is ok to cast away the constness of cin.*/

    // Downcast to int
    int Mint = (int) M;

    p_veryend = LTFAT_FFTW(plan_many_dft)(1, &Mint, N * W,
                                          (LTFAT_COMPLEX *)cin, NULL,
                                          1, Mint,
                                          cwork, NULL,
                                          1, Mint,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);

    /* -------- Execute initial IFFT ------------------------ */
    LTFAT_FFTW(execute)(p_veryend);


    if (ptype)
    {
        LTFAT_NAME_COMPLEX(dgtphaseunlockhelper)(cwork, L, W, a, M, cwork);

       /* for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt n = 0; n < N; n++)
            {
                LTFAT_COMPLEX* cworktmp = cwork + w * N * M + n * M;
                LTFAT_NAME_COMPLEX(circshift)(cworktmp, cworktmp, M, a * n);
            }
        }
        */
    }


    /* -------- Main loop ----------------------------------- */

    const ltfatInt ld4c = M * N;

    /* Leading dimensions of cf */
    const ltfatInt ld3b = q * q * W;

    /* Leading dimensions of the 4dim array. */
    const ltfatInt ld2ff = p * q * W;

    for (ltfatInt r = 0; r < c; r++)
    {


        LTFAT_COMPLEX *cfp = cf;

        for (ltfatInt w = 0; w < W; w++)
        {
            /* Complete inverse fac of coefficients */
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt u = 0; u < q; u++)
                {
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        const ltfatInt rem = r + l * c +
                                             positiverem(u + s * q - l * h_a, N)
                                             * M + w * ld4c;
                        cbuf[s] = cwork[rem];
                    }

                    /* Do inverse fft of length d */
                    LTFAT_FFTW(execute)(p_after);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        cfp[s * ld3b] =  cbuf[s];
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
        for (ltfatInt s = 0; s < d; s++)
        {

            const LTFAT_COMPLEX *gbase = gf + (r + s * c) * p * q;
            LTFAT_COMPLEX       *fbase = ff + s * p * q * W;
            const LTFAT_COMPLEX *cbase = (const LTFAT_COMPLEX *)cf + s * q * q * W;

            for (ltfatInt nm = 0; nm < q * W; nm++)
            {
                for (ltfatInt km = 0; km < p; km++)
                {

                    fbase[km + nm * p] = 0.0;
                    for (ltfatInt mm = 0; mm < q; mm++)
                    {
                        fbase[km + nm * p] += gbase[km + mm * p] * cbase[mm + nm * q];
                    }
                    /* Scale because of FFTWs normalization. */
                    fbase[km + nm * p] = fbase[km + nm * p] * scalconst;
                }
            }
        }




        /* ----------- compute inverse signal factorization ---------- */


        LTFAT_COMPLEX *ffp = ff;
        LTFAT_COMPLEX *fp  = f + r;

        for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt k = 0; k < p; k++)
                {
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        cbuf[s] = ffp[s * ld2ff];
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        const ltfatInt rem = positiverem(k * M + s * p * M -
                                                         l * h_a * a, L);
                        fp[rem] = cbuf[s];
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
        fp -= L * W;

        /* ----- Main loop ends here ------------- */
    }

    /* -----------  Clean up ----------------- */

    LTFAT_FFTW(destroy_plan)(p_veryend);
    LTFAT_FFTW(destroy_plan)(p_after);
    LTFAT_FFTW(destroy_plan)(p_before);

    LTFAT_SAFEFREEALL(cwork, ff, cf, cbuf);
}



LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf,
                         const ltfatInt L, const ltfatInt W,
                         const ltfatInt a, const ltfatInt M,
                         const dgt_phasetype ptype, LTFAT_REAL *f)
{

    /*  --------- initial declarations -------------- */

    ltfatInt h_a, h_m;

    LTFAT_FFTW(plan) p_before, p_after, p_veryend;
    LTFAT_COMPLEX *ff, *cf, *cbuf;
    LTFAT_REAL *cwork, *sbuf;

    /* This is a floor operation. */
    const ltfatInt M2 = M / 2 + 1;

    /*  ----------- calculation of parameters and plans -------- */

    const ltfatInt b = L / M;
    const ltfatInt N = L / a;

    const ltfatInt c = gcd(a, M, &h_a, &h_m);
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = b / p;

    /* This is a floor operation. */
    const ltfatInt d2 = d / 2 + 1;

    h_a = -h_a;

    ff    = ltfat_malloc(d2 * p * q * W * sizeof * ff);
    cf    = ltfat_malloc(d2 * q * q * W * sizeof * cf);
    cwork = ltfat_malloc(M * N * W * sizeof * cwork);
    cbuf  = ltfat_malloc(d2 * sizeof * cbuf);
    sbuf  = ltfat_malloc(d * sizeof * sbuf);

    /* Scaling constant needed because of FFTWs normalization. */
    const LTFAT_REAL scalconst = 1.0 / ((LTFAT_REAL)d * sqrt((LTFAT_REAL)M));

    /* Create plans. In-place. */
    p_before = LTFAT_FFTW(plan_dft_c2r_1d)(d, cbuf, sbuf, FFTW_ESTIMATE);

    p_after  = LTFAT_FFTW(plan_dft_r2c_1d)(d, sbuf, cbuf, FFTW_ESTIMATE);

    /* Create plan. Copy data so we do not overwrite input. Therefore
       it is ok to cast away the constness of cin. This transform
       destroys its input by default, but the extra flag should prevent
       this. */

    // Downcast to int
    int Mint = (int) M;
    p_veryend = LTFAT_FFTW(plan_many_dft_c2r)(1, &Mint, N * W,
                (LTFAT_COMPLEX *)cin, NULL,
                1, M2,
                cwork, NULL,
                1, M,
                FFTW_ESTIMATE + FFTW_PRESERVE_INPUT);


    /* -------- Execute initial IFFT ------------------------ */
    LTFAT_FFTW(execute)(p_veryend);

    if (ptype)
    {

        LTFAT_NAME_REAL(dgtphaseunlockhelper)(cwork, L, W, a, M, cwork);
        /*for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt n = 0; n < N; n++)
            {
                LTFAT_REAL* cworktmp = cwork + w * N * M + n * N;
                LTFAT_NAME_REAL(circshift)(cworktmp, cworktmp, N, a * n);
            }
        }
        */
    }


    const ltfatInt ld4c = M * N;

    /* Leading dimensions of cf */
    const ltfatInt ld3b = q * q * W;

    /* Leading dimensions of the 4dim array. */
    const ltfatInt ld2ff = p * q * W;

    /* -------- Main loop ----------------------------------- */
    for (ltfatInt r = 0; r < c; r++)
    {

        /* -------- compute coefficient factorization ----------- */

        LTFAT_COMPLEX *cfp = cf;

        for (ltfatInt w = 0; w < W; w++)
        {
            /* Complete inverse fac of coefficients */
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt u = 0; u < q; u++)
                {
                    for (ltfatInt s = 0; s < d; s++)
                    {
                        sbuf[s] = cwork[r + l * c +
                                        positiverem(u + s * q - l * h_a, N) *
                                        M + w * ld4c];
                    }

                    /* Do inverse fft of length d */
                    LTFAT_FFTW(execute)(p_after);

                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        cfp[s * ld3b] = cbuf[s];
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
            const LTFAT_COMPLEX *gbase = gf + (r + s * c) * p * q;
            LTFAT_COMPLEX       *fbase = ff + s * p * q * W;
            const LTFAT_COMPLEX *cbase = (const LTFAT_COMPLEX *)cf + s * q * q * W;

            for (ltfatInt nm = 0; nm < q * W; nm++)
            {
                for (ltfatInt km = 0; km < p; km++)
                {
                    fbase[km + nm * p] = 0.0;
                    for (ltfatInt mm = 0; mm < q; mm++)
                    {
                        fbase[km + nm * p] += gbase[km + mm * p] * cbase[mm + nm * q];
                    }
                    /* Scale because of FFTWs normalization. */
                    fbase[km + nm * p] = fbase[km + nm * p] * scalconst;
                }
            }
        }


        /* ----------- compute inverse signal factorization ---------- */

        LTFAT_COMPLEX *ffp = ff;
        LTFAT_REAL    *fp  = f + r;

        for (ltfatInt w = 0; w < W; w++)
        {
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt k = 0; k < p; k++)
                {
                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        cbuf[s] = ffp[s * ld2ff];
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        fp[positiverem(k * M + s * p * M - l * h_a * a, L)] = sbuf[s];
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
        fp -= L * W;


        /*  ------- main loop ends -------- */
    }

    /* -----------  Clean up ----------------- */

    LTFAT_FFTW(destroy_plan)(p_veryend);
    LTFAT_FFTW(destroy_plan)(p_after);
    LTFAT_FFTW(destroy_plan)(p_before);

    LTFAT_SAFEFREEALL(cwork, ff, cf, cbuf, sbuf);
}
