#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(fftcircshift)( const LTFAT_COMPLEX* in, int L, double shift,
                          LTFAT_COMPLEX* out)
{
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
}


LTFAT_EXTERN void
LTFAT_NAME(fftfftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out)
{
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        LTFAT_NAME(fftcircshift)(in, L, domod.quot, out);
    }
    else
    {
        out[0] = in[0];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < L; ii += 2)
            out[ii] = -1.0 * in[ii];
    }
}

LTFAT_EXTERN void
LTFAT_NAME(fftifftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out)
{
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        LTFAT_NAME(fftcircshift)(in, L, -domod.quot, out);
    else
        LTFAT_NAME(fftfftshift)(in, L, out);
}

LTFAT_EXTERN void
LTFAT_NAME(fftrealcircshift)( const LTFAT_COMPLEX* in, int L, double shift,
                              LTFAT_COMPLEX* out)
{
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
}

LTFAT_EXTERN void
LTFAT_NAME(fftrealfftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out)
{
    const div_t domod = div(L, 2);

    if (domod.rem)
    {
        // There is no Nyquist sample, modulation is by (L-1/L)*pi
        LTFAT_NAME(fftrealcircshift)(in, L, domod.quot, out);
    }
    else
    {
        out[0] = in[0];

        // There is Nyquist sample. Modulation is exactly by pi i.e. no complex
        // multiplication
        for (int ii = 1; ii < domod.quot + 1; ii += 2)
            out[ii] = -1.0 * in[ii];
    }
}

LTFAT_EXTERN void
LTFAT_NAME(fftrealifftshift)(const LTFAT_COMPLEX* in, int L, LTFAT_COMPLEX* out)
{
    const div_t domod = div(L, 2);

    if (domod.rem)
        // There is no Nyquist sample, modulation is by -(L-1/L)*pi
        LTFAT_NAME(fftrealcircshift)(in, L, -domod.quot, out);
    else
        LTFAT_NAME(fftrealfftshift)(in, L, out);
}
