#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(circshift)(const LTFAT_TYPE in[], const ltfatInt L,
                      const ltfatInt shift, LTFAT_TYPE out[])
{
    // Fix shift
    ltfatInt p = (L - shift) % L;

    if (p < 0) p += L;

    if (in == out)
    {
        ltfatInt m, count, i, j;

        // Circshift inplace is magic!
        for (m = 0, count = 0; count != L; m++)
        {
            LTFAT_TYPE t = in[m];

            for (i = m, j = m + p; j != m;
                 i = j, j = j + p < L ? j + p : j + p - L, count++)
                out[i] = out[j];

            out[i] = t; count++;
        }
    }
    else
    {
        // Circshit out of place is boring ...
        memcpy(out, in + p, (L - p)*sizeof * out);
        memcpy(out + L - p, in, p * sizeof * out);
    }
}


// in might be equal to out
LTFAT_EXTERN void
LTFAT_NAME(fftshift)(const LTFAT_TYPE* in, ltfatInt L, LTFAT_TYPE* out)
{
    LTFAT_NAME(circshift)(in, L, (L / 2), out);
}

// in might be equal to out
LTFAT_EXTERN void
LTFAT_NAME(ifftshift)(const LTFAT_TYPE* in, ltfatInt L, LTFAT_TYPE* out)
{
    LTFAT_NAME(circshift)(in, L, -(L / 2), out);
}


LTFAT_EXTERN void
LTFAT_NAME(reverse_array)(LTFAT_TYPE* in, LTFAT_TYPE* out,
                          const ltfatInt L)
{

    if (in == out)
    {
        LTFAT_TYPE tmpVar = (LTFAT_TYPE) 0.0;

        for (ltfatInt ii = 0; ii < L / 2; ii++)
        {
            tmpVar = in[L - 1 - ii];
            in[L - 1 - ii] = in[ii];
            in[ii] = tmpVar;
        }
    }
    else
    {
        for (ltfatInt ii = 0; ii < L; ii++)
        {
            out[ii] = in[L - 1 - ii];
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(conjugate_array)(LTFAT_TYPE* in, LTFAT_TYPE* out,
                            const ltfatInt L)
{
#ifdef LTFAT_COMPLEXTYPE

    for (ltfatInt ii = 0; ii < L; ii++)
    {
        out[ii] = LTFAT_COMPLEXH(conj)(in[ii]);
    }

#else

    if (in == out)
    {
        return;
    }
    else
    {
        memcpy(out, in, L * sizeof(LTFAT_TYPE));
    }

#endif

}

LTFAT_EXTERN void
LTFAT_NAME(periodize_array)(LTFAT_TYPE* in, const ltfatInt Lin,
                            LTFAT_TYPE* out, const ltfatInt Lout)
{
    /* Do nothing if there is no place where to put periodized samples */
    if ( Lout <= Lin )
    {
        if ( in != out )
        {
            memcpy(out, in, Lout * sizeof * in);
        }

        return;
    }

    ltfatInt periods = floor( Lout / ((double)Lin) );
    ltfatInt lastL = Lout - periods * Lin;
    ltfatInt startPer = in == out ? 1 : 0;

    for (ltfatInt ii = startPer; ii < periods; ii++)
    {
        memcpy(out + ii * Lin, in, Lin * sizeof * in);
    }

    memcpy(out + periods * Lin, in, lastL * sizeof * in);
}

LTFAT_EXTERN void
LTFAT_NAME(fold_array)(const LTFAT_TYPE* in, const ltfatInt Lin,
                       const ltfatInt Lfold, LTFAT_TYPE* out)
{
    if (Lfold > Lin)
    {
        // Just copy
        if (in != out)
        {
            memcpy(out, in, Lin * sizeof * out);
            memset(out + Lin, 0, (Lfold - Lin)*sizeof * out);
        }

        return;
    }

    ltfatInt startIdx = 0;

    if (in != out)
        memset(out, 0, Lfold * sizeof * out);
    else
        startIdx = Lfold;

    for (ltfatInt ii = startIdx; ii < Lin; ii++)
        out[ii % Lfold] += in[ii];
}


LTFAT_EXTERN void
LTFAT_NAME(array2complex)(LTFAT_TYPE* in, LTFAT_COMPLEX* out,
                          const ltfatInt L)
{
#ifdef LTFAT_COMPLEXTYPE

    if (in == (LTFAT_TYPE*)out)
    {
        return;
    }
    else
    {
        memcpy(out, in, L * sizeof(LTFAT_COMPLEX));
    }

#else

    if (in == (LTFAT_TYPE*)out)
    {
        // This should produce an error
    }
    else
    {
        LTFAT_REAL (*outTmp)[2] = (LTFAT_REAL(*)[2])  out;

        for (ltfatInt ii = 0; ii < L; ii++)
        {
            outTmp[ii][0] = in[ii];
            outTmp[ii][1] = (LTFAT_TYPE) 0.0;
        }
    }

#endif
}

LTFAT_EXTERN void
LTFAT_NAME(dgtphaselockhelper)(LTFAT_TYPE* cin, const ltfatInt L,
                               const ltfatInt W, const ltfatInt a,
                               const ltfatInt M, LTFAT_TYPE* cout)
{
    ltfatInt N = L / a;

    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt n = 0; n < N; n++)
        {
            ltfatInt offset = w * N * M + n * M;
            LTFAT_TYPE* cintmp = cin + offset;
            LTFAT_TYPE* couttmp = cout + offset;
            // circshift takes care of possible inplace operation
            LTFAT_NAME(circshift)(cintmp, M, -a * n, couttmp);
        }

    }

}

LTFAT_EXTERN void
LTFAT_NAME(dgtphaseunlockhelper)(LTFAT_TYPE* cin, const ltfatInt L,
                                 const ltfatInt W, const ltfatInt a,
                                 const ltfatInt M, LTFAT_TYPE* cout)
{
    ltfatInt N = L / a;

    for (ltfatInt w = 0; w < W; w++)
    {
        for (ltfatInt n = 0; n < N; n++)
        {
            ltfatInt offset = w * N * M + n * M;
            LTFAT_TYPE* cintmp = cin + offset;
            LTFAT_TYPE* couttmp = cout + offset;
            // circshift takes care of possible inplace operation
            LTFAT_NAME(circshift)(cintmp, M, a * n, couttmp);
        }

    }

}

LTFAT_EXTERN void
LTFAT_NAME(findmaxinarray)(const LTFAT_TYPE* in, const ltfatInt L,
                           LTFAT_TYPE* max, ltfatInt* idx)
{
    *max = in[0];
    *idx = 0;

    for (ltfatInt ii = 1; ii < L; ++ii)
    {
#ifdef LTFAT_COMPLEXTYPE

        if (LTFAT_COMPLEXH(cabs)(in[ii]) > LTFAT_COMPLEXH(cabs)(*max) )
#else
        if (in[ii] > *max)
#endif
        {
            *max = in[ii];
            *idx = ii;
        }
    }
}

LTFAT_EXTERN int
LTFAT_NAME(findmaxinarraywrtmask)(const LTFAT_TYPE* in, const int* mask,
                                  const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx)
{
    int found = 0;
    *max = -1e99;
    *idx = 0;

    for (ltfatInt ii = 0; ii < L; ++ii)
    {
#ifdef LTFAT_COMPLEXTYPE

        if (!mask[ii] && LTFAT_COMPLEXH(cabs)(in[ii]) > LTFAT_COMPLEXH(cabs)(*max))
#else
        if (!mask[ii] && in[ii] > *max)
#endif
        {
            *max = in[ii];
            *idx = ii;
            found = 1;
        }
    }

    return found;
}

/* This routine changes an FIR window to a LONG window.
 *
 * Input parameters:
 *  f     : Real valued input array.
 *  Lfir  : Length of input array
 *  Llong  : Length of output array
 *  h     : Output array
 */
LTFAT_EXTERN void
LTFAT_NAME(fir2long)(const LTFAT_TYPE* f, const ltfatInt Lfir,
                     const ltfatInt Llong, LTFAT_TYPE* h)
{
    const div_t domod = div(Lfir, 2);

    /* ---- In the odd case, the additional element is kept in the first half. ---*/

    for (ltfatInt ii = 0; ii < domod.quot + domod.rem; ii++)
    {
        h[ii] = f[ii];
    }

    for (ltfatInt ii = domod.quot + domod.rem; ii < Llong - domod.quot; ii++)
    {
        h[ii] = 0.0;
    }

    const ltfatInt ss = Llong - Lfir;

    for (ltfatInt ii = Lfir - 1; ii >= domod.quot + domod.rem; ii--)
    {
        h[ii + ss] = f[ii];
    }

}

/* This routine changes a LONG window to an FIR window.
 *
 * Input parameters:
 *  f     : Real valued input array.
 *  Llong  : Length of input array
 *  Lfir  : Length of output array
 *  h     : Output array
 */
LTFAT_EXTERN void
LTFAT_NAME(long2fir)(const LTFAT_TYPE* f, const ltfatInt Llong,
                     const ltfatInt Lfir, LTFAT_TYPE* h)
{
    const div_t domod = div(Lfir, 2);

    /* ---- In the odd case, the additional element is kept in the first half. ---*/

    for (ltfatInt ii = 0; ii < domod.quot + domod.rem; ii++)
    {
        h[ii] = f[ii];
    }

    const ltfatInt ss = Llong - Lfir;

    for (ltfatInt ii = domod.quot + domod.rem; ii < Lfir; ii++)
    {
        h[ii] = f[ii + ss];
    }

}

LTFAT_EXTERN void
LTFAT_NAME(normalize)(const LTFAT_TYPE* in, const ltfatInt L,
                      ltfat_normalize_t flag, LTFAT_TYPE* out)
{
    LTFAT_REAL normfac = 1.0;

    switch (flag)
    {
    case LTFAT_NORMALIZE_ENERGY:
    {
        normfac = 0.0;

        for (ltfatInt ii = 0; ii < L; ii++)
            normfac += in[ii] * in[ii];

        normfac = sqrt(normfac);
    }

    case LTFAT_NORMALIZE_NULL:
    default:
        normfac = 1.0;
    };

    for (ltfatInt ii = 0; ii < L; ii++)
        out[ii] = normfac * in[ii];
}
