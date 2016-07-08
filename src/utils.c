#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

LTFAT_EXTERN int
LTFAT_NAME(fftcircshift)( const LTFAT_COMPLEX* in, const ltfatInt L,
                          const double shift, LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");

    if (shift != 0)
    {
        LTFAT_COMPLEX phasefact = -I * 2.0 * M_PI * shift / L;

        out[0] = in[0];

        for (int ii = 1; ii < L; ii++)
            out[ii] = cexp(ii * phasefact) * in[ii];
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
LTFAT_NAME(fftfftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                        LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        CHECKSTATUS(LTFAT_NAME(fftcircshift)(in, L, domod.quot, out),
                    "fftcircshift failed");
    }
    else
    {
        CHECKNULL(in); CHECKNULL(out);
        CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");

        out[0] = in[0];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < L; ii += 2)
            out[ii] = -1.0 * in[ii];
    }
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(fftifftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                         LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        CHECKSTATUS(LTFAT_NAME(fftcircshift)(in, L, -domod.quot, out),
                    "fftcircshift failed");
    else
        CHECKSTATUS(LTFAT_NAME(fftfftshift)(in, L, out),
                    "fftcircshift failed");
error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(fftrealcircshift)( const LTFAT_COMPLEX* in, const ltfatInt L,
                              const double shift, LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");

    const div_t domod = div(L, 2);

    if (shift != 0)
    {
        LTFAT_COMPLEX phasefact = -I * 2.0 * M_PI * shift / L;

        out[0] = in[0];

        for (int ii = 1; ii < domod.quot + 1; ii++)
            out[ii] = cexp(ii * phasefact) * in[ii];
    }
    else
    {
        if (in != out)
            memcpy(out, in, domod.quot * sizeof * out);
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(fftrealfftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                            LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        CHECKSTATUS(LTFAT_NAME(fftrealcircshift)(in, L, domod.quot, out),
                    "fftrealcircshift failed");
    }
    else
    {
        CHECKNULL(in); CHECKNULL(out);
        CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");

        out[0] = in[0];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < domod.quot + 1; ii += 2)
            out[ii] = -1.0 * in[ii];
    }

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(fftrealifftshift)(const LTFAT_COMPLEX* in, const ltfatInt L,
                             LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        CHECKSTATUS(LTFAT_NAME(fftrealcircshift)(in, L, -domod.quot, out),
                    "fftrealcircshift failed");
    else
        CHECKSTATUS(LTFAT_NAME(fftrealfftshift)(in, L, out),
                    "fftrealcircshift");

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(real2complex_array)(const LTFAT_REAL* in, const ltfatInt L,
                               LTFAT_COMPLEX* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");

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
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");

    const LTFAT_REAL (*inTmp)[2] = (const LTFAT_REAL(*)[2]) in;

    for (ltfatInt ii = 0; ii < L; ii++)
        out[ii] = inTmp[ii][0];

error:
    return status;
}
