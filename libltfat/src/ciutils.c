#include "ltfat.h"
#include "ltfat_types.h"
#include "ciutils.h"

LTFAT_EXTERN
void LTFAT_NAME(circshift)(LTFAT_TYPE* in, LTFAT_TYPE* out, const ltfatInt L,
                           const ltfatInt shift)
{
    ltfatInt shiftMod = shift % L;

    if (in == out)
    {

        if (1)
        {
            LTFAT_TYPE* inTmp = (LTFAT_TYPE*)ltfat_malloc(L * sizeof(LTFAT_TYPE));
            memcpy(inTmp, in, L * sizeof(LTFAT_TYPE));
            LTFAT_NAME(circshift)(inTmp, out, L, shift);
            ltfat_free(inTmp);
        }
        else
        {
            ltfatInt m, count, ii, jj;

            for (m = 0, count = 0; count != L; m++)
            {
                LTFAT_TYPE t = in[m];

                for (ii = m, jj = m + shiftMod;
                     jj != m;
                     ii = jj, jj = jj + shiftMod < L ? jj + shiftMod : jj + shiftMod - L, count++)
                {
                    in[ii] = in[jj];
                }

                in[ii] = t;
                count++;
            }
        }



        return;
    }



    if (shiftMod < 0)
    {
        memcpy(out, in - shiftMod, (L + shiftMod)*sizeof * out);
        memcpy(out + (L + shiftMod), in, -shiftMod * sizeof * out);
    }
    else if (shiftMod > 0)
    {
        memcpy(out + shiftMod, in, (L - shiftMod)*sizeof * out);
        memcpy(out, in + L - shiftMod, shiftMod * sizeof * out);
    }
    else
    {
        memcpy(out, in, L * sizeof * out);
    }
}

LTFAT_EXTERN
void LTFAT_NAME(reverse_array)(LTFAT_TYPE* in, LTFAT_TYPE* out,
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

LTFAT_EXTERN
void LTFAT_NAME(conjugate_array)(LTFAT_TYPE* in, LTFAT_TYPE* out,
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

LTFAT_EXTERN
void LTFAT_NAME(periodize_array)(LTFAT_TYPE* in, const ltfatInt Lin,
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

LTFAT_EXTERN
void LTFAT_NAME(array2complex)(LTFAT_TYPE* in, LTFAT_COMPLEX* out,
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
            LTFAT_NAME(circshift)(cintmp, couttmp, M, -a * n);
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
            LTFAT_NAME(circshift)(cintmp, couttmp, M, a * n);
        }

    }

}

LTFAT_EXTERN
void LTFAT_NAME(findmaxinarray)(const LTFAT_TYPE* in, const ltfatInt L,
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

LTFAT_EXTERN
int LTFAT_NAME(findmaxinarraywrtmask)(const LTFAT_TYPE* in, const int* mask,
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

