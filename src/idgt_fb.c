#include "ltfat.h"
#include "ltfat_types.h"

#define THE_SUM { \
     /* Copy to c-buffer and ifft it */ \
     for (ltfatInt m=0; m<M; m++) \
     { \
        cbuf[m]=cin[m+n*M+w*M*N]; \
     } \
     LTFAT_FFTW(execute)(p_small); \
\
    ltfatInt premarg = ptype?glh:-n*a+glh;  \
     const ltfatInt rem = positiverem(premarg, M); \
     for (ltfatInt ii=0; ii<gl/M; ii++) \
     { \
        for (ltfatInt m=0; m<rem; m++) \
        { \
           ff[m+ii*M]=cbuf[M-rem+m]*gw[m+ii*M]; \
        } \
        for (ltfatInt m=0; m<M-rem; m++) \
        { \
           ff[m+ii*M+rem] = cbuf[m]*gw[m+rem+ii*M]; \
        } \
     } \
}

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
                    const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                    const ltfatInt a, const ltfatInt M,
                    const dgt_phasetype ptype, LTFAT_COMPLEX *f)

{
    /*  --------- initial declarations -------------- */

    const ltfatInt N = L / a;

    ltfatInt ep, sp;

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));
    LTFAT_COMPLEX *fw;

    LTFAT_COMPLEX *cbuf = ltfat_malloc(M * sizeof * cbuf);

    /* Create plan. In-place. */
    LTFAT_FFTW(plan) p_small = LTFAT_FFTW(plan_dft_1d)(M, cbuf, cbuf,
                               FFTW_BACKWARD, FFTW_MEASURE);

    /* % The fftshift actually makes some things easier. */
    LTFAT_COMPLEX *gw  = ltfat_malloc(gl * sizeof * gw);
    for (ltfatInt l = 0; l < glh; l++)
    {
        gw[l] = g[l + (gl - glh)];
    }
    for (ltfatInt l = glh; l < gl; l++)
    {
        gw[l] = g[l - glh];
    }

    LTFAT_COMPLEX *ff  = ltfat_malloc(gl * sizeof * ff);

    for (ltfatInt w = 0; w < W; w++)
    {
        fw = f + w * L;
        for (ltfatInt l = 0; l < L; l++)
        {
            fw[l] = (LTFAT_COMPLEX) 0.0;
        }
        /* ----- Handle the first boundary using periodic boundary conditions. --- */
        for (ltfatInt n = 0; n < glh_d_a; n++)
        {

            THE_SUM;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* % Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
            {
                fw[sp + ii] += ff[ii];
            }
            for (ltfatInt ii = 0; ii < ep + 1; ii++)
            {
                fw[ii] += ff[L - sp + ii];
            }
        }


        /* ----- Handle the middle case. --------------------- */
        for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
        {

            THE_SUM;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < ep - sp + 1; ii++)
            {
                fw[ii + sp] += ff[ii];
            }
        }

        /* Handle the last boundary using periodic boundary conditions. */
        for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
        {

            THE_SUM;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
            {
                fw[sp + ii] += ff[ii];
            }
            for (ltfatInt ii = 0; ii < ep + 1; ii++)
            {
                fw[ii] += ff[L - sp + ii];
            }

        }
    }

    LTFAT_SAFEFREEALL(cbuf, ff, gw);
    LTFAT_FFTW(destroy_plan)(p_small);

}

/* ------------------- IDGT_FB_R --------------------- */


#define THE_SUM_R { \
     \
     for (ltfatInt m=0; m<M; m++) \
     { \
        cbuf[m]=cin[m+n*M+w*M*N]; \
     } \
     LTFAT_FFTW(execute)(p_small); \
\
     const ltfatInt rem = positiverem(-n*a+glh, M); \
     for (ltfatInt ii=0; ii<gl/M; ii++) \
     { \
        for (ltfatInt m=0; m<rem; m++) \
        { \
           ff[m+ii*M] = cbuf[M-rem+m]*gw[m+ii*M]; \
        } \
        for (ltfatInt m=0; m<M-rem; m++) \
        { \
           ff[m+ii*M+rem] = cbuf[m]*gw[m+rem+ii*M]; \
        } \
     } \
}

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb_r)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                      const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                      const ltfatInt a, const ltfatInt M,
                      LTFAT_COMPLEX *f)

{
    /*  --------- initial declarations -------------- */

    const ltfatInt N = L / a;

    ltfatInt ep, sp;

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));
    LTFAT_COMPLEX *fw;

    LTFAT_COMPLEX *cbuf = (LTFAT_COMPLEX*)ltfat_malloc(M * sizeof(LTFAT_COMPLEX));

    /* Create plan. In-place. */
    LTFAT_FFTW(plan) p_small = LTFAT_FFTW(plan_dft_1d)(M, cbuf, cbuf,
                               FFTW_BACKWARD, FFTW_MEASURE);

    /* % The fftshift actually makes some things easier. */
    LTFAT_REAL *gw  = (LTFAT_REAL*)ltfat_malloc(gl * sizeof(LTFAT_REAL));
    for (ltfatInt l = 0; l < glh; l++)
    {
        gw[l] = g[l + (gl - glh)];
    }
    for (ltfatInt l = glh; l < gl; l++)
    {
        gw[l] = g[l - glh];
    }

    LTFAT_COMPLEX *ff  = (LTFAT_COMPLEX*)ltfat_malloc(gl * sizeof(LTFAT_COMPLEX));

    for (ltfatInt w = 0; w < W; w++)
    {
        fw = f + w * L;
        for (ltfatInt l = 0; l < L; l++)
        {
            fw[l] = (LTFAT_COMPLEX)0.0;
        }
        /* ----- Handle the first boundary using periodic boundary conditions. --- */
        for (ltfatInt n = 0; n < glh_d_a; n++)
        {

            THE_SUM_R;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* % Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
            {
                fw[sp + ii] += ff[ii];
            }
            for (ltfatInt ii = 0; ii < ep + 1; ii++)
            {
                fw[ii] += ff[L - sp + ii];
            }
        }


        /* ----- Handle the middle case. --------------------- */
        for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
        {

            THE_SUM_R;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < ep - sp + 1; ii++)
            {
                fw[ii + sp] += ff[ii];
            }
        }

        /* Handle the last boundary using periodic boundary conditions. */
        for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
        {

            THE_SUM_R;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
            {
                fw[sp + ii] += ff[ii];
            }
            for (ltfatInt ii = 0; ii < ep + 1; ii++)
            {
                fw[ii] += ff[L - sp + ii];
            }

        }
    }

    LTFAT_SAFEFREEALL(cbuf, ff, gw);
    LTFAT_FFTW(destroy_plan)(p_small);

}




/* ------------------- IDGTREAL ---------------------- */

#define THE_SUM_REAL { \
     for (ltfatInt m=0; m<M2; m++) \
     { \
        cbuf[m]=cin[m+n*M2+w*M2*N]; \
     } \
     LTFAT_FFTW(execute)(p_small); \
    \
    ltfatInt premarg = ptype?glh:-n*a+glh;  \
     const ltfatInt rem = positiverem(premarg, M); \
     for (ltfatInt ii=0; ii<gl/M; ii++) \
     { \
        for (ltfatInt m=0; m<rem; m++) \
        { \
           ff[m+ii*M] = crbuf[M-rem+m]*gw[m+ii*M]; \
        } \
        for (ltfatInt m=0; m<M-rem; m++) \
        { \
           ff[m+ii*M+rem] = crbuf[m]*gw[m+rem+ii*M]; \
        } \
     } \
}

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                        const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, LTFAT_REAL *f)

{
    /*  --------- initial declarations -------------- */

    const ltfatInt N = L / a;

    /* This is a floor operation. */
    const ltfatInt M2 = M / 2 + 1;

    ltfatInt ep, sp;

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));
    LTFAT_REAL *fw;

    LTFAT_COMPLEX *cbuf  = ltfat_malloc(M2 * sizeof * cbuf);
    LTFAT_REAL    *crbuf = ltfat_malloc( M * sizeof * crbuf);

    /* Create plan. In-place. */
    LTFAT_FFTW(plan) p_small = LTFAT_FFTW(plan_dft_c2r_1d)(M, cbuf, crbuf, FFTW_MEASURE);

    /* % The fftshift actually makes some things easier. */
    LTFAT_REAL *gw  = ltfat_malloc(gl * sizeof * gw);
    for (ltfatInt l = 0; l < glh; l++)
    {
        gw[l] = g[l + (gl - glh)];
    }
    for (ltfatInt l = glh; l < gl; l++)
    {
        gw[l] = g[l - glh];
    }

    LTFAT_REAL *ff  = ltfat_malloc(gl * sizeof * ff);

    for (ltfatInt w = 0; w < W; w++)
    {
        fw = f + w * L;
        for (ltfatInt l = 0; l < L; l++)
        {
            fw[l] = 0.0;
        }
        /* ----- Handle the first boundary using periodic boundary conditions. --- */
        for (ltfatInt n = 0; n < glh_d_a; n++)
        {

            THE_SUM_REAL;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* % Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
            {
                fw[sp + ii] += ff[ii];
            }
            for (ltfatInt ii = 0; ii < ep + 1; ii++)
            {
                fw[ii] += ff[L - sp + ii];
            }
        }


        /* ----- Handle the middle case. --------------------- */
        for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
        {

            THE_SUM_REAL;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < ep - sp + 1; ii++)
            {
                fw[ii + sp] += ff[ii];
            }
        }

        /* Handle the last boundary using periodic boundary conditions. */
        for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
        {

            THE_SUM_REAL;

            sp = positiverem(n * a - glh, L);
            ep = positiverem(n * a - glh + gl - 1, L);

            /* Add the ff vector to f at position sp. */
            for (ltfatInt ii = 0; ii < L - sp; ii++)
            {
                fw[sp + ii] += ff[ii];
            }
            for (ltfatInt ii = 0; ii < ep + 1; ii++)
            {
                fw[ii] += ff[L - sp + ii];
            }

        }
    }

    LTFAT_SAFEFREEALL(cbuf, crbuf, ff, gw);
    LTFAT_FFTW(destroy_plan)(p_small);

}
