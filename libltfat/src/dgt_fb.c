#include "ltfat.h"
#include "ltfat_types.h"

#define CH(name) LTFAT_COMPLEXH(name)

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
for (ltfatInt m=0;m<M;m++) \
{ \
   ltfatInt premarg = plan.ptype?-glh:n*a-glh; \
   rem = 2*positiverem(m+(premarg), M); \
   sbuf[rem]=0.0; \
   sbuf[rem+1]=0.0; \
   fbd=fw+2*m; \
   for (ltfatInt k=0;k<gl/M;k++) \
   { \
      sbuf[rem]  += fbd[0]; \
      sbuf[rem+1]+= fbd[1]; \
      fbd+=2*M; \
   } \
} \
\
 LTFAT_FFTW(execute)(plan.p_small);         \
\
coefsum=(LTFAT_REAL*)cout+2*(n*M+w*M*N); \
for (ltfatInt m=0;m<M;m++) \
{ \
   coefsum[2*m] = sbuf[2*m]; \
   coefsum[2*m+1] = sbuf[2*m+1]; \
}}

#define THE_SUM_REAL { \
ltfatInt premarg = plan.ptype?-glh:n*a-glh; \
for (ltfatInt m=0;m<M;m++) \
{ \
   rem = positiverem(m+(premarg), M); \
   sbuf[rem]=0.0; \
   fbd=fw+m; \
   for (ltfatInt k=0;k<gl/M;k++) \
   { \
     sbuf[rem]+=(*fbd);         \
      fbd+=M; \
   } \
} \
\
 LTFAT_FFTW(execute)(plan.p_small);         \
\
coefsum=(LTFAT_REAL*)cout+2*(n*M2+w*M2*N); \
for (ltfatInt m=0;m<M2;m++) \
{ \
   coefsum[2*m]   = CH(creal)(cbuf[m]); \
   coefsum[2*m+1] = CH(cimag)(cbuf[m]); \
}}


LTFAT_EXTERN LTFAT_NAME(dgt_fb_plan)
LTFAT_NAME(dgt_fb_init)(const LTFAT_COMPLEX *g,
                        const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, unsigned flags)
{
    LTFAT_NAME(dgt_fb_plan) plan;

    plan.a = a;
    plan.M = M;
    plan.gl = gl;
    plan.ptype = ptype;

    plan.gw  = (LTFAT_COMPLEX*)ltfat_malloc(gl * sizeof(LTFAT_COMPLEX));

    plan.fw  = (LTFAT_REAL*)ltfat_malloc(2 * gl * sizeof(LTFAT_REAL));

    plan.sbuf = (LTFAT_REAL*)ltfat_malloc(M * sizeof(LTFAT_COMPLEX));

    plan.p_small = LTFAT_FFTW(plan_dft_1d)(M, (LTFAT_COMPLEX*)plan.sbuf,
                                           (LTFAT_COMPLEX*)plan.sbuf,
                                           FFTW_FORWARD, flags);

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;


    /* Do the fftshift of g to place the center in the middle and
     * conjugate it.
     */

    for (ltfatInt l = 0; l < glh; l++)
    {
        plan.gw[l] = CH(conj)(g[l + (gl - glh)]);
    }
    for (ltfatInt l = glh; l < gl; l++)
    {
        plan.gw[l] = CH(conj)(g[l - glh]);
    }

    return (plan);
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_done)(LTFAT_NAME(dgt_fb_plan) plan)
{
    LTFAT_SAFEFREEALL(plan.sbuf, plan.gw, plan.fw);
    LTFAT_FFTW(destroy_plan)(plan.p_small);
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_execute)(LTFAT_NAME(dgt_fb_plan) plan, const LTFAT_COMPLEX *f,
                           const ltfatInt L, const ltfatInt W,  LTFAT_COMPLEX *cout)
{
    /*  --------- initial declarations -------------- */

    const ltfatInt a = plan.a;
    const ltfatInt M = plan.M;
    const ltfatInt N = L / a;

    const ltfatInt gl = plan.gl;
    LTFAT_REAL *sbuf = plan.sbuf;
    LTFAT_REAL *fw = plan.fw;

    /* This is a floor operation. */
    const ltfatInt glh = plan.gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));

    ltfatInt rem;

    LTFAT_COMPLEX *gb;
    LTFAT_REAL *coefsum, *fbd;


    /*  ---------- main body ----------- */

    /*----- Handle the first boundary using periodic boundary conditions.*/
    for (ltfatInt n = 0; n < glh_d_a; n++)
    {
        gb = plan.gw;
        for (ltfatInt w = 0; w < W; w++)
        {

            fbd = (LTFAT_REAL*)f + 2 * (L - (glh - n * a) + L * w);
            for (ltfatInt l = 0; l < glh - n * a; l++)
            {
                fw[2 * l]  = fbd[2 * l] * CH(creal)(gb[l]) - fbd[2 * l + 1] * CH(cimag)(gb[l]);
                fw[2 * l + 1] = fbd[2 * l + 1] * CH(creal)(gb[l]) + fbd[2 * l] * CH(cimag)(gb[l]);
            }
            fbd = (LTFAT_REAL*)f - 2 * (glh - n * a) + 2 * L * w;
            for (ltfatInt l = glh - n * a; l < gl; l++)
            {
                fw[2 * l]  = fbd[2 * l] * CH(creal)(gb[l]) - fbd[2 * l + 1] * CH(cimag)(gb[l]);
                fw[2 * l + 1] = fbd[2 * l + 1] * CH(creal)(gb[l]) + fbd[2 * l] * CH(cimag)(gb[l]);
            }

            THE_SUM

        }
    }

    /* ----- Handle the middle case. --------------------- */
    for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
    {
        gb = plan.gw;
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = (LTFAT_REAL*)f + 2 * (n * a - glh + L * w);
            for (ltfatInt l = 0; l < gl; l++)
            {
                fw[2 * l]  = fbd[2 * l] * CH(creal)(gb[l]) - fbd[2 * l + 1] * CH(cimag)(gb[l]);
                fw[2 * l + 1] = fbd[2 * l + 1] * CH(creal)(gb[l]) + fbd[2 * l] * CH(cimag)(gb[l]);
            }

            THE_SUM
        }

    }

    /* Handle the last boundary using periodic boundary conditions. */
    for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
    {
        gb = plan.gw;
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = (LTFAT_REAL*)f + 2 * (n * a - glh + L * w);
            for (ltfatInt l = 0; l < L - n * a + glh; l++)
            {
                fw[2 * l]  = fbd[2 * l] * CH(creal)(gb[l]) - fbd[2 * l + 1] * CH(cimag)(gb[l]);
                fw[2 * l + 1] = fbd[2 * l + 1] * CH(creal)(gb[l]) + fbd[2 * l] * CH(cimag)(gb[l]);
            }
            fbd = (LTFAT_REAL*)f - 2 * (L - n * a + glh) + 2 * L * w;
            for (ltfatInt l = L - n * a + glh; l < gl; l++)
            {
                fw[2 * l]  = fbd[2 * l] * CH(creal)(gb[l]) - fbd[2 * l + 1] * CH(cimag)(gb[l]);
                fw[2 * l + 1] = fbd[2 * l + 1] * CH(creal)(gb[l]) + fbd[2 * l] * CH(cimag)(gb[l]);
            }

            THE_SUM
        }
    }

}

/* See the comments on the macro THE_SUM. This macro uses real valued
 * inputs, but produces complex valued output and uses a regular FFT.
 */
#define THE_SUM_R {for (ltfatInt m=0;m<M;m++) \
{ \
   ltfatInt premarg = ptype?-glh:n*a-glh; \
   rem = 2*positiverem(m+(premarg), M); \
   sbuf[rem]=0.0; \
   sbuf[rem+1]=0.0; \
   fbd=fw+m; \
   for (ltfatInt k=0;k<gl/M;k++) \
   { \
      sbuf[rem]  += (*fbd); \
      fbd+=M; \
   } \
}   \
\
    LTFAT_FFTW(execute)(p_small);       \
\
coefsum=(LTFAT_REAL*)cout+2*(n*M+r*M*N+w*M*N*R); \
for (ltfatInt m=0;m<M;m++) \
{ \
   coefsum[2*m]   = sbuf[2*m]; \
   coefsum[2*m+1] = sbuf[2*m+1]; \
}}


LTFAT_EXTERN void
LTFAT_NAME(dgt_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                   const ltfatInt L, const ltfatInt gl,
                   const ltfatInt W, const ltfatInt a, const ltfatInt M,
                  const dgt_phasetype ptype, LTFAT_COMPLEX *cout)
{
    /*  --------- initial declarations -------------- */

    ltfatInt r, rem;

    LTFAT_REAL *gw;

    LTFAT_FFTW(plan) p_small;

    LTFAT_REAL *gb;
    LTFAT_REAL *sbuf, *coefsum, *fw;

    const LTFAT_REAL *fbd;


    const ltfatInt R = 1;

    /*  ----------- calculation of parameters and plans -------- */

    const ltfatInt N = L / a;

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));

    gw   = (LTFAT_REAL*)ltfat_malloc(gl * R * sizeof(LTFAT_REAL));
    fw   = (LTFAT_REAL*)ltfat_malloc(gl * sizeof(LTFAT_REAL));
    sbuf = (LTFAT_REAL*)ltfat_malloc(2 * M * sizeof(LTFAT_REAL));

    /* Create plan. In-place. */
    p_small = LTFAT_FFTW(plan_dft_1d)(M, (LTFAT_COMPLEX*)sbuf,
                                      (LTFAT_COMPLEX*)sbuf,
                                      FFTW_FORWARD, FFTW_MEASURE);

    /*  ---------- main body ----------- */

    /* Do the fftshift of g to place the center in the middle and
     * conjugate it.
     */

    for (r = 0; r < R; r++)
    {
        for (ltfatInt l = 0; l < glh; l++)
        {
            gw[l + gl * r] = g[l + (gl - glh) + gl * r];
        }
        for (ltfatInt l = glh; l < gl; l++)
        {
            gw[l + gl * r] = g[l - glh + gl * r];
        }
    }

    /*----- Handle the first boundary using periodic boundary conditions.*/
    for (ltfatInt n = 0; n < glh_d_a; n++)
    {
        for (ltfatInt r = 0; r < R; r++)
        {
            gb = gw + r * gl;
            for (ltfatInt w = 0; w < W; w++)
            {

                fbd = f + L - (glh - n * a) + L * w;
                for (ltfatInt l = 0; l < glh - n * a; l++)
                {
                    fw[l]  = fbd[l] * gb[l];
                }
                fbd = f - (glh - n * a) + L * w;
                for (ltfatInt l = glh - n * a; l < gl; l++)
                {
                    fw[l]  = fbd[l] * gb[l];
                }

                THE_SUM_R

            }

        }
    }

    /* ----- Handle the middle case. --------------------- */
    for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
    {

        for (ltfatInt r = 0; r < R; r++)
        {
            gb = gw + r * gl;
            for (ltfatInt w = 0; w < W; w++)
            {
                fbd = f + (n * a - glh + L * w);
                for (ltfatInt l = 0; l < gl; l++)
                {
                    fw[l]  = fbd[l] * gb[l];
                }

                THE_SUM_R
            }

        }

    }

    /* Handle the last boundary using periodic boundary conditions. */
    for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
    {
        for (ltfatInt r = 0; r < R; r++)
        {
            gb = gw + r * gl;
            for (ltfatInt w = 0; w < W; w++)
            {
                fbd = f + (n * a - glh + L * w);
                for (ltfatInt l = 0; l < L - n * a + glh; l++)
                {
                    fw[l]  = fbd[l] * gb[l];
                }
                fbd = f - (L - n * a + glh) + L * w;
                for (ltfatInt l = L - n * a + glh; l < gl; l++)
                {
                    fw[l]  = fbd[l] * gb[l];
                }

                THE_SUM_R
            }
        }
    }

    /* -----------  Clean up ----------------- */
    LTFAT_SAFEFREEALL(sbuf, gw, fw);
    LTFAT_FFTW(destroy_plan)(p_small);
}








LTFAT_EXTERN LTFAT_NAME(dgtreal_fb_plan)
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL *g,
                            const ltfatInt gl, const ltfatInt a,
                            const ltfatInt M, const dgt_phasetype ptype,
                            unsigned flags)
{
    LTFAT_NAME(dgtreal_fb_plan) plan;

    plan.a = a;
    plan.M = M;
    const ltfatInt M2 = M / 2 + 1;
    plan.gl = gl;
    plan.ptype = ptype;


    plan.gw   = ltfat_malloc(gl * sizeof * plan.gw);

    plan.fw   = ltfat_malloc(gl * sizeof * plan.fw);

    plan.sbuf = ltfat_malloc(M * sizeof * plan.sbuf);

    plan.cbuf = ltfat_malloc(M2 * sizeof * plan.cbuf);

    plan.p_small = LTFAT_FFTW(plan_dft_r2c_1d)(M, plan.sbuf, plan.cbuf, flags);

    /* This is a floor operation. */
    const ltfatInt glh = gl / 2;


    /* Do the fftshift of g to place the center in the middle and
     * conjugate it.
     */

    for (ltfatInt l = 0; l < glh; l++)
    {
        plan.gw[l] = g[l + (gl - glh)];
    }
    for (ltfatInt l = glh; l < gl; l++)
    {
        plan.gw[l] = g[l - glh];
    }

    return (plan);
}

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan) plan)
{
    LTFAT_SAFEFREEALL(plan.sbuf, plan.cbuf, plan.gw, plan.fw);
    LTFAT_FFTW(destroy_plan)(plan.p_small);
}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_execute)(const LTFAT_NAME(dgtreal_fb_plan) plan,
                               const LTFAT_REAL *f,
                               const ltfatInt L, const ltfatInt W,
                               LTFAT_COMPLEX *cout)
{
    /*  --------- initial declarations -------------- */

    const ltfatInt a = plan.a;
    const ltfatInt M = plan.M;
    const ltfatInt N = L / a;

    const ltfatInt gl = plan.gl;
    LTFAT_REAL    *sbuf = plan.sbuf;
    LTFAT_COMPLEX *cbuf = plan.cbuf;
    LTFAT_REAL *fw = plan.fw;

    /* These are floor operations. */
    const ltfatInt glh = plan.gl / 2;
    const ltfatInt M2 = M / 2 + 1;

    /* This is a ceil operation. */
    const ltfatInt glh_d_a = (ltfatInt)ceil((glh * 1.0) / (a));

    ltfatInt rem;

    LTFAT_REAL *gb;
    LTFAT_REAL *coefsum;

    const LTFAT_REAL *fbd;

    /*  ---------- main body ----------- */

    /*----- Handle the first boundary using periodic boundary conditions.*/
    for (ltfatInt n = 0; n < glh_d_a; n++)
    {
        gb = plan.gw;
        for (ltfatInt w = 0; w < W; w++)
        {

            fbd = f + L - (glh - n * a) + L * w;
            for (ltfatInt l = 0; l < glh - n * a; l++)
            {
                fw[l]  = fbd[l] * gb[l];
            }
            fbd = f - (glh - n * a) + L * w;
            for (ltfatInt l = glh - n * a; l < gl; l++)
            {
                fw[l]  = fbd[l] * gb[l];
            }

            THE_SUM_REAL

        }
    }

    /* ----- Handle the middle case. --------------------- */
    for (ltfatInt n = glh_d_a; n < (L - (gl + 1) / 2) / a + 1; n++)
    {

        gb = plan.gw;
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = f + (n * a - glh + L * w);
            for (ltfatInt l = 0; l < gl; l++)
            {
                fw[l]  = fbd[l] * gb[l];
            }

            THE_SUM_REAL
        }

    }

    /* Handle the last boundary using periodic boundary conditions. */
    for (ltfatInt n = (L - (gl + 1) / 2) / a + 1; n < N; n++)
    {
        gb = plan.gw;
        for (ltfatInt w = 0; w < W; w++)
        {
            fbd = f + (n * a - glh + L * w);
            for (ltfatInt l = 0; l < L - n * a + glh; l++)
            {
                fw[l]  = fbd[l] * gb[l];
            }
            fbd = f - (L - n * a + glh) + L * w;
            for (ltfatInt l = L - n * a + glh; l < gl; l++)
            {
                fw[l]  = fbd[l] * gb[l];
            }

            THE_SUM_REAL
        }
    }

}

#undef THE_SUM
#undef CH
