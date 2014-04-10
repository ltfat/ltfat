#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN LTFAT_NAME(dgt_ola_plan)
LTFAT_NAME(dgt_ola_init)(const LTFAT_COMPLEX *g, const ltfatInt gl,
                         const ltfatInt W, const ltfatInt a, const ltfatInt M,
                         const ltfatInt bl, const dgt_phasetype ptype,
                         unsigned flags)
{

    LTFAT_NAME(dgt_ola_plan) plan;

    plan.bl = bl;
    plan.gl = gl;
    plan.W  = W;

    const ltfatInt Lext    = bl + gl;
    const ltfatInt Nblocke = Lext / a;

    plan.buf  = ltfat_malloc(Lext * W * sizeof * plan.buf);
    plan.gext = ltfat_malloc(Lext * sizeof * plan.gext);
    plan.cbuf = ltfat_malloc(M * Nblocke * W * sizeof * plan.cbuf);

    LTFAT_NAME(fir2long_c)(g, gl, Lext, plan.gext);

    /* Zero the last part of the buffer, it will always be zero. */
    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt jj = bl; jj < Lext; jj++)
        {
            plan.buf[jj + w * Lext] = (LTFAT_COMPLEX) 0.0;
        }
    }

    plan.plan =
        LTFAT_NAME(dgt_long_init)((const LTFAT_COMPLEX*)plan.buf,
                                  (const LTFAT_COMPLEX*)plan.gext,
                                  Lext, W, a, M,
                                  plan.cbuf, ptype, flags);

    return (plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_execute)(const LTFAT_NAME(dgt_ola_plan) plan,
                            const LTFAT_COMPLEX *f, const ltfatInt L,
                            LTFAT_COMPLEX *cout)

{
    const ltfatInt bl      = plan.bl;
    const ltfatInt gl      = plan.gl;
    const ltfatInt a       = plan.plan.a;
    const ltfatInt M       = plan.plan.M;
    const ltfatInt N       = L / a;
    const ltfatInt Lext    = bl + gl;
    const ltfatInt Nb      = L / bl;
    const ltfatInt b2      = gl / a / 2;
    const ltfatInt Nblock  = bl / a;
    const ltfatInt Nblocke = Lext / a;
    const ltfatInt W       = plan.W;


    /* Zero the output array, as we will be adding to it */
    for (ltfatInt ii = 0; ii < M * N * W; ii++)
    {
        cout[ii] = (LTFAT_COMPLEX) 0.0;
    }

    for (ltfatInt ii = 0; ii < Nb; ii++)
    {
        ltfatInt s_ii;

        /* Copy to working buffer. */
        for (ltfatInt w = 0; w < W; w++)
        {
            memcpy(plan.buf + Lext * w, f + ii * bl + w * L, sizeof(LTFAT_COMPLEX)*bl);
        }

        /* Execute the short DGT */
        LTFAT_NAME(dgt_long_execute)(plan.plan);

        /* Place the results */
        for (ltfatInt w = 0; w < W; w++)
        {
            LTFAT_COMPLEX *cout_p;
            LTFAT_COMPLEX *cbuf_p;

            /* Place large block */
            cout_p = cout + ii * M * Nblock + w * M * N ;
            cbuf_p = plan.cbuf +             w * M * Nblocke;
            for (ltfatInt m = 0; m < M; m++)
            {
                for (ltfatInt n = 0; n < Nblock; n++)
                {
                    cout_p[m + n * M] += cbuf_p[m + n * M];
                }
            }

            /* Small block + */
            s_ii = positiverem(ii + 1, Nb);
            cout_p = cout + s_ii * M * Nblock + w * M * N ;
            cbuf_p = plan.cbuf +      M * Nblock + w * M * Nblocke;
            for (ltfatInt m = 0; m < M; m++)
            {
                for (ltfatInt n = 0; n < b2; n++)
                {
                    cout_p[m + n * M] += cbuf_p[m + n * M];
                }
            }


            /* Small block - */
            s_ii = positiverem(ii - 1, Nb) + 1;
            cout_p = cout + M * (s_ii * Nblock - b2) + w * M * N ;
            cbuf_p = plan.cbuf + M * (Nblock + b2)     + w * M * Nblocke;
            for (ltfatInt m = 0; m < M; m++)
            {
                for (ltfatInt n = 0; n < b2; n++)
                {
                    cout_p[m + n * M] += cbuf_p[m + n * M];
                }
            }

        }

    }


}

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_done)(LTFAT_NAME(dgt_ola_plan) plan)
{
    LTFAT_NAME(dgt_long_done)(plan.plan);
    LTFAT_SAFEFREEALL(plan.cbuf, plan.gext, plan.buf);
}



LTFAT_EXTERN LTFAT_NAME(dgtreal_ola_plan)
LTFAT_NAME(dgtreal_ola_init)(const LTFAT_REAL *g, const ltfatInt gl,
                             const ltfatInt W, const ltfatInt a, const ltfatInt M,
                             const ltfatInt bl, const dgt_phasetype ptype,
                             unsigned flags)
{

    LTFAT_NAME(dgtreal_ola_plan) plan;

    plan.bl = bl;
    plan.gl = gl;
    plan.W  = W;
    const ltfatInt M2 = M / 2 + 1;

    const ltfatInt Lext    = bl + gl;
    const ltfatInt Nblocke = Lext / a;

    plan.buf  = (LTFAT_REAL*) ltfat_malloc(Lext * W * sizeof(LTFAT_REAL));
    plan.gext = (LTFAT_REAL*) ltfat_malloc(Lext * sizeof(LTFAT_REAL));
    plan.cbuf = (LTFAT_COMPLEX*) ltfat_malloc(M2 * Nblocke * W * sizeof(LTFAT_COMPLEX));

    LTFAT_NAME(fir2long_r)(g, gl, Lext, plan.gext);

    /* Zero the last part of the buffer, it will always be zero. */
    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt jj = bl; jj < Lext; jj++)
        {
            plan.buf[jj + w * Lext] = 0.0;
        }
    }

    plan.plan =
        LTFAT_NAME(dgtreal_long_init)((const LTFAT_REAL*)plan.buf,
                                      (const LTFAT_REAL*)plan.gext,
                                      Lext, W, a, M, plan.cbuf,
                                      ptype, flags);

    return (plan);

}







LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_execute)(const LTFAT_NAME(dgtreal_ola_plan) plan,
                                const LTFAT_REAL *f, const ltfatInt L,
                                LTFAT_COMPLEX *cout)

{
    const ltfatInt bl      = plan.bl;
    const ltfatInt gl      = plan.gl;
    const ltfatInt a       = plan.plan.a;
    const ltfatInt M       = plan.plan.M;
    const ltfatInt N       = L / a;
    const ltfatInt Lext    = bl + gl;
    const ltfatInt Nb      = L / bl;
    const ltfatInt b2      = gl / a / 2;
    const ltfatInt Nblock  = bl / a;
    const ltfatInt Nblocke = Lext / a;
    const ltfatInt W       = plan.W;
    const ltfatInt M2      = M / 2 + 1;

    /* Zero the output array, as we will be adding to it */
    for (ltfatInt ii = 0; ii < M2 * N * W; ii++)
    {
        cout[ii] = (LTFAT_COMPLEX) 0.0;
    }


    for (ltfatInt ii = 0; ii < Nb; ii++)
    {
        ltfatInt s_ii;

        /* Copy to working buffer. */
        for (ltfatInt w = 0; w < W; w++)
        {
            memcpy(plan.buf + Lext * w, f + ii * bl + w * L, sizeof(LTFAT_REAL)*bl);
        }

        /* Execute the short DGTREAL */
        LTFAT_NAME(dgtreal_long_execute)(plan.plan);

        /* Place the results */
        for (ltfatInt w = 0; w < W; w++)
        {
            LTFAT_COMPLEX *cout_p;
            LTFAT_COMPLEX *cbuf_p;

            /* Place large block */
            cout_p = cout + ii * M2 * Nblock + w * M2 * N ;
            cbuf_p = plan.cbuf +             w * M2 * Nblocke;
            for (ltfatInt m = 0; m < M2; m++)
            {
                for (ltfatInt n = 0; n < Nblock; n++)
                {
                    cout_p[m + n * M2] += cbuf_p[m + n * M2];
                }
            }

            /* Small block + */
            s_ii = positiverem(ii + 1, Nb);
            cout_p = cout + s_ii * M2 * Nblock + w * M2 * N ;
            cbuf_p = plan.cbuf +      M2 * Nblock + w * M2 * Nblocke;
            for (ltfatInt m = 0; m < M2; m++)
            {
                for (ltfatInt n = 0; n < b2; n++)
                {
                    cout_p[m + n * M2] += cbuf_p[m + n * M2];
                }
            }


            /* Small block - */
            s_ii = positiverem(ii - 1, Nb) + 1;
            cout_p = cout + M2 * (s_ii * Nblock - b2) + w * M2 * N ;
            cbuf_p = plan.cbuf + M2 * (     Nblock + b2) + w * M2 * Nblocke;
            for (ltfatInt m = 0; m < M2; m++)
            {
                for (ltfatInt n = 0; n < b2; n++)
                {
                    cout_p[m + n * M2] += cbuf_p[m + n * M2];
                }
            }

        }

    }
}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_done)(LTFAT_NAME(dgtreal_ola_plan) plan)
{
    LTFAT_NAME(dgtreal_long_done)(plan.plan);
    LTFAT_SAFEFREEALL(plan.cbuf, plan.gext, plan.buf);
}
