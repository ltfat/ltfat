#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(fftcircshift)( const LTFAT_COMPLEX* in, const ltfatInt L,
                                  const double shift, LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    double shiftLoc = remainder(shift, L);

    if (shiftLoc != 0)
    {
        const div_t domod = div(L, 2);
        LTFAT_COMPLEX phasefact = -I * 2.0 * M_PI * shiftLoc / L;

        out[0] = in[0];

        for (int ii = 1; ii < domod.quot + 1; ii++)
            out[ii] = LTFAT_COMPLEXH(cexp)(ii * phasefact) * in[ii];

        for (int ii = 1; ii < domod.quot + domod.rem; ii++)
            out[L-ii] = LTFAT_COMPLEXH(cexp)(-ii * phasefact) * in[L-ii];

    }
    else
    {
        if (in != out)
            memcpy(out, in, L * sizeof * out);
    }
error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(fftfftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                                LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        status = LTFAT_NAME_COMPLEX(fftcircshift)(in, L, domod.quot, out);
    }
    else
    {
        CHECKNULL(in); CHECKNULL(out);
        CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

        if (in != out)
            for (int ii = 0; ii < L; ii += 2)
                out[ii] = in[ii];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < L; ii += 2)
            out[ii] = -1.0 * in[ii];
    }
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(fftifftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                                 LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        status = LTFAT_NAME_COMPLEX(fftcircshift)(in, L, -domod.quot, out);
    else
        status = LTFAT_NAME_COMPLEX(fftfftshift)(in, L, out);
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(fftrealcircshift)( const LTFAT_COMPLEX* in, const ltfatInt L,
                                      const double shift, LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    const div_t domod = div(L, 2);

    double shiftLoc = remainder(shift, L);

    if (shiftLoc != 0)
    {
        LTFAT_COMPLEX phasefact = -I * 2.0 * M_PI * shiftLoc / L;

        out[0] = in[0];

        for (int ii = 1; ii < domod.quot + 1; ii++)
            out[ii] = LTFAT_COMPLEXH(cexp)(ii * phasefact) * in[ii];
    }
    else
    {
        if (in != out)
            memcpy(out, in, (domod.quot + 1) * sizeof * out);
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(fftrealfftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                                    LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        status = LTFAT_NAME_COMPLEX(fftrealcircshift)(in, L, domod.quot, out);
    }
    else
    {
        CHECKNULL(in); CHECKNULL(out);
        CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

        if (in != out)
            for (int ii = 0; ii < domod.quot + 1; ii += 2)
                out[ii] = in[ii];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < domod.quot + 1; ii += 2)
            out[ii] = -1.0 * in[ii];
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(fftrealifftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                                     LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        status = LTFAT_NAME_COMPLEX(fftrealcircshift)(in, L, -domod.quot, out);
    else
        status = LTFAT_NAME_COMPLEX(fftrealfftshift)(in, L, out);

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(real2complex_array)(const LTFAT_REAL* in, const ltfatInt L,
                               LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    LTFAT_REAL (*outTmp)[2] = (LTFAT_REAL(*)[2])  out;

    if (in == (LTFAT_REAL*)out)
    {
        // Go from the back to avoid overwriting input
        for (ltfatInt ii = L - 1; ii >= 0; ii--)
        {
            outTmp[ii][0] = in[ii];
            outTmp[ii][1] = 0.0;
        }
    }
    else
    {
        for (ltfatInt ii = 0; ii < L; ii++)
        {
            outTmp[ii][0] = in[ii];
            outTmp[ii][1] = 0.0;
        }
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(complex2real_array)(const LTFAT_COMPLEX* in, const ltfatInt L,
                               LTFAT_REAL* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    const LTFAT_REAL (*inTmp)[2] = (const LTFAT_REAL(*)[2]) in;

    for (ltfatInt ii = 0; ii < L; ii++)
        out[ii] = inTmp[ii][0];

error:
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(dgt_phaselock)(const LTFAT_COMPLEX* cFreqinv,
                                  const ltfatInt L, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                                  LTFAT_COMPLEX* cTimeinv)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cFreqinv); CHECKNULL(cTimeinv);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    ltfatInt N = L / a;

    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt n = 0; n < N; n++)
        {
            const LTFAT_COMPLEX* inCol = cFreqinv + n * M + w * M * N;
            LTFAT_COMPLEX* outCol = cTimeinv + n * M + w * M * N;
            LTFAT_NAME_COMPLEX(fftcircshift)( inCol, M, -n * a, outCol);
        }
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(dgtreal_phaselock)(const LTFAT_COMPLEX* cFreqinv,
                                      const ltfatInt L, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                                      LTFAT_COMPLEX* cTimeinv)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cFreqinv); CHECKNULL(cTimeinv);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    ltfatInt N = L / a;
    ltfatInt M2 = M / 2 + 1;

    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt n = 0; n < N; n++)
        {
            const LTFAT_COMPLEX* inCol = cFreqinv + n * M2 + w * M2 * N;
            LTFAT_COMPLEX* outCol = cTimeinv + n * M2 + w * M2 * N;
            LTFAT_NAME_COMPLEX(fftrealcircshift)( inCol, M, -n * a, outCol);
        }
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(dgt_phaseunlock)(const LTFAT_COMPLEX* cTimeinv,
                                    const ltfatInt L, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                                    LTFAT_COMPLEX* cFreqinv)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cFreqinv); CHECKNULL(cTimeinv);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    ltfatInt N = L / a;

    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt n = 0; n < N; n++)
        {
            const LTFAT_COMPLEX* inCol = cTimeinv + n * M + w * M * N;
            LTFAT_COMPLEX* outCol = cFreqinv + n * M + w * M * N;
            LTFAT_NAME_COMPLEX(fftcircshift)( inCol, M, n * a, outCol);
        }
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME_COMPLEX(dgtreal_phaseunlock)(const LTFAT_COMPLEX* cTimeinv,
                                        const ltfatInt L, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                                        LTFAT_COMPLEX* cFreqinv)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(cFreqinv); CHECKNULL(cTimeinv);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    ltfatInt N = L / a;
    ltfatInt M2 = M / 2 + 1;

    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt n = 0; n < N; n++)
        {
            const LTFAT_COMPLEX* inCol = cTimeinv + n * M2 + w * M2 * N;
            LTFAT_COMPLEX* outCol = cFreqinv + n * M2 + w * M2 * N;
            LTFAT_NAME_COMPLEX(fftrealcircshift)( inCol, M, n * a, outCol);
        }
    }

error:
    return status;
}
