#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

// in might be equal to out
LTFAT_API int
LTFAT_NAME(circshift)(const LTFAT_TYPE in[], ltfat_int L,
                      ltfat_int shift, LTFAT_TYPE out[])
{
    ltfat_int p;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    // Fix shift
    p = (L - shift) % L;

    if (p < 0) p += L;

    if (in == out)
    {
        if (p) // Do nothing if no shift is needed
        {
            ltfat_int m, count, i, j;

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
    }
    else
    {
        // Still ok if p==0
        memcpy(out, in + p, (L - p)*sizeof * out);
        memcpy(out + L - p, in, p * sizeof * out);
    }

error:
    return status;
}


// in might be equal to out
LTFAT_API int
LTFAT_NAME(fftshift)(const LTFAT_TYPE* in, ltfat_int L, LTFAT_TYPE* out)
{
    return LTFAT_NAME(circshift)(in, L, (L / 2), out);
}

// in might be equal to out
LTFAT_API int
LTFAT_NAME(ifftshift)(const LTFAT_TYPE* in, ltfat_int L, LTFAT_TYPE* out)
{
    return LTFAT_NAME(circshift)(in, L, -(L / 2), out);
}


LTFAT_API int
LTFAT_NAME(reverse_array)(const LTFAT_TYPE* in, ltfat_int L,
                          LTFAT_TYPE* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    if (in == out)
    {
        LTFAT_TYPE tmpVar = (LTFAT_TYPE) 0.0;

        for (ltfat_int ii = 0; ii < L / 2; ii++)
        {
            tmpVar = out[L - 1 - ii];
            out[L - 1 - ii] = out[ii];
            out[ii] = tmpVar;
        }
    }
    else
    {
        for (ltfat_int ii = 0; ii < L; ii++)
            out[ii] = in[L - 1 - ii];
    }

error:
    return status;
}

LTFAT_API int
LTFAT_NAME(conjugate_array)(const LTFAT_TYPE* in, ltfat_int L,
                            LTFAT_TYPE* out)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

#ifdef LTFAT_COMPLEXTYPE
    for (ltfat_int ii = 0; ii < L; ii++)
        out[ii] = conj(in[ii]); // type-generic macro conj
#else
    if (in != out)
        memcpy(out, in, L * sizeof * out);
#endif

error:
    return status;

}

LTFAT_API int
LTFAT_NAME(periodize_array)(const LTFAT_TYPE* in, ltfat_int Lin,
                            ltfat_int Lout, LTFAT_TYPE* out )
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, Lin > 0, "Lin must be positive");
    CHECK(LTFATERR_BADSIZE, Lout > 0, "Lout must be positive");

    /* Do nothing if there is no place where to put periodized samples */
    if ( Lout <= Lin )
    {
        if ( in != out )
            memcpy(out, in, Lout * sizeof * in);
    }
    else
    {
        ltfat_int periods =  Lout / Lin;
        ltfat_int lastL = Lout - periods * Lin;
        ltfat_int startPer = in == out ? 1 : 0;

        for (ltfat_int ii = startPer; ii < periods; ii++)
            memcpy(out + ii * Lin, in, Lin * sizeof * in);

        memcpy(out + periods * Lin, in, lastL * sizeof * in);
    }
error:
    return status;
}

/*
  *
 * If offset is not zero, the function performs:
 * fold_array(in,Lin,0,Lfold,out);
 * circshift(out,Lfold,offset,out);
 *
 * or equivalently
 *
 * circshift(in,Lin,offset,in);
 * fold_array(in,Lin,0,Lfold,out);
 *
 * without the intermediate step.
 * */
LTFAT_API int
LTFAT_NAME(fold_array)(const LTFAT_TYPE* in, ltfat_int Lin,
                       ltfat_int offset,
                       ltfat_int Lfold, LTFAT_TYPE* out)
{
    ltfat_int startIdx;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, Lin > 0, "Lin must be positive");
    CHECK(LTFATERR_BADSIZE, Lfold > 0, "Lfold must be positive");

    // Sanitize offset.
    startIdx = ltfat_positiverem(offset, Lfold);

    // Clear output, we will use it as an accumulator
    if (in != out)
        memset(out, 0, Lfold * sizeof * out);
    else if (Lfold > Lin)
        memset(out + Lin, 0, (Lfold - Lin)*sizeof * out);

    if (!startIdx)
    {
        // Common code for no offset
        ltfat_int startAt = in == out ? Lfold : 0;

        for (ltfat_int ii = startAt; ii < Lin;)
            for (ltfat_int kk = 0; ii < Lin && kk < Lfold; ii++, kk++)
                out[kk] += in[ii];

    }
    else
    {
        if (in == out)
        {
            // We cannot avoid the (slow) inplace circshift anyway.
            // Lets do it after the folding
            LTFAT_NAME(fold_array)(in, Lin, 0, Lfold, out);
            LTFAT_NAME(circshift)(in, Lfold, startIdx, out);
        }
        else
        {
            // We avoid the inplace circshift by effectivelly
            // doing circshift of all blocks
            for (ltfat_int ii = 0; ii < Lin;)
            {
                ltfat_int kk = startIdx;
                for (; kk < Lfold && ii < Lin; ii++, kk++)
                    out[kk] += in[ii];

                for (kk = 0; kk < startIdx && ii < Lin; ii++, kk++)
                    out[kk] += in[ii];
            }
        }
    }
error:
    return status;
}


LTFAT_API int
LTFAT_NAME(ensurecomplex_array)(const LTFAT_TYPE* in,  ltfat_int L,
                                LTFAT_COMPLEX* out)
{
#ifdef LTFAT_COMPLEXTYPE
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    if (in != (LTFAT_TYPE*)out)
        memcpy(out, in, L * sizeof * out);

error:
    return status;
#else
    return LTFAT_NAME_REAL(real2complex_array)(in, L, out);
#endif
}


LTFAT_API void
LTFAT_NAME(dgtphaselockhelper)(LTFAT_TYPE* cin, ltfat_int L,
                               ltfat_int W, ltfat_int a,
                               ltfat_int M, LTFAT_TYPE* cout)
{
    ltfat_int N = L / a;

    for (ltfat_int w = 0; w < W; w++)
    {
        for (ltfat_int n = 0; n < N; n++)
        {
            ltfat_int offset = w * N * M + n * M;
            LTFAT_TYPE* cintmp = cin + offset;
            LTFAT_TYPE* couttmp = cout + offset;
            LTFAT_NAME(circshift)(cintmp, M, -a * n, couttmp);
        }

    }

}

LTFAT_API void
LTFAT_NAME(dgtphaseunlockhelper)(LTFAT_TYPE* cin, ltfat_int L,
                                 ltfat_int W, ltfat_int a,
                                 ltfat_int M, LTFAT_TYPE* cout)
{
    ltfat_int N = L / a;

    for (ltfat_int w = 0; w < W; w++)
    {
        for (ltfat_int n = 0; n < N; n++)
        {
            ltfat_int offset = w * N * M + n * M;
            LTFAT_TYPE* cintmp = cin + offset;
            LTFAT_TYPE* couttmp = cout + offset;
            LTFAT_NAME(circshift)(cintmp, M, a * n, couttmp);
        }

    }

}

LTFAT_API void
LTFAT_NAME(findmaxinarray)(const LTFAT_TYPE* in, ltfat_int L,
                           LTFAT_TYPE* max, ltfat_int* idx)
{
    *max = in[0];
    *idx = 0;

    for (ltfat_int ii = 1; ii < L; ++ii)
    {
#ifdef LTFAT_COMPLEXTYPE

        if ( ltfat_abs(in[ii]) > ltfat_abs(*max) )
#else
        if (in[ii] > *max)
#endif
        {
            *max = in[ii];
            *idx = ii;
        }
    }
}

LTFAT_API int
LTFAT_NAME(findmaxinarraywrtmask)(const LTFAT_TYPE* in, const int* mask,
                                  ltfat_int L, LTFAT_TYPE* max, ltfat_int* idx)
{
    int found = 0;
    *max = (LTFAT_REAL) -1e99;
    *idx = 0;

    for (ltfat_int ii = 0; ii < L; ++ii)
    {
#ifdef LTFAT_COMPLEXTYPE
        if (!mask[ii] && ltfat_abs(in[ii]) > ltfat_abs(*max))
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

LTFAT_API int
LTFAT_NAME(fir2long)(const LTFAT_TYPE* in, ltfat_int Lfir,
                     ltfat_int Llong, LTFAT_TYPE* out)
{
    ltfat_div_t domod;
    ltfat_int ss;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, Llong > 0, "Llong must be positive");
    CHECK(LTFATERR_BADSIZE, Lfir > 0, "Lfir must be positive");
    CHECK(LTFATERR_BADREQSIZE, Lfir <= Llong, "Lfir <= Llong does not hold");

    domod = ltfat_idiv(Lfir, 2);

    /* ---- In the odd case, the additional element is kept in the first half. ---*/

    // Copy first half
    if (in != out)
        memcpy(out, in, (domod.quot + domod.rem)*sizeof * out);

    ss = Llong - Lfir;
    // Copy second half from the back
    for (ltfat_int ii = Lfir - 1; ii >= domod.quot + domod.rem; ii--)
        out[ii + ss] = in[ii];

    // Zero out the middle
    for (ltfat_int ii = domod.quot + domod.rem; ii < Llong - domod.quot ; ii++)
        out[ii] = 0.0;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(long2fir)(const LTFAT_TYPE* in, ltfat_int Llong,
                     ltfat_int Lfir, LTFAT_TYPE* out)
{
    ltfat_div_t domod;
    ltfat_int ss;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, Llong > 0, "Llong must be positive");
    CHECK(LTFATERR_BADSIZE, Lfir > 0, "Lfir must be positive");
    CHECK(LTFATERR_BADREQSIZE, Lfir <= Llong, "Lfir <= Llong does not hold");

    domod = ltfat_idiv(Lfir, 2);

    /* ---- In the odd case, the additional element is kept in the first half. ---*/

    ss = Llong - Lfir;

    if (in != out)
        memcpy(out, in, (domod.quot + domod.rem)*sizeof * out);

    for (ltfat_int ii = domod.quot + domod.rem; ii < Lfir; ii++)
        out[ii] = in[ii + ss];

error:
    return status;
}

LTFAT_API int
LTFAT_NAME(normalize)(const LTFAT_TYPE* in, ltfat_int L,
                      ltfat_normalize_t flag, LTFAT_TYPE* out)
{
    LTFAT_REAL normfac;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive");

    normfac = 1.0;

    switch (flag)
    {
    case LTFAT_NORMALIZE_ENERGY:
    {
        normfac = 0.0;

        for (ltfat_int ii = 0; ii < L; ii++)
        {
#ifdef LTFAT_COMPLEXTYPE
            LTFAT_REAL inAbs = ltfat_abs(in[ii]);
#else
            LTFAT_REAL inAbs = in[ii]; // We dont need abs here
#endif
            normfac += inAbs * inAbs;
        }

        normfac = sqrt(normfac);
        break;
    }
    case LTFAT_NORMALIZE_AREA:
    {
        normfac = 0.0;

        for (ltfat_int ii = 0; ii < L; ii++)
        {
            LTFAT_REAL inAbs = ltfat_abs(in[ii]);
            normfac += inAbs;
        }


        break;
    }
    case LTFAT_NORMALIZE_PEAK:
    {
        normfac = 0.0;

        for (ltfat_int ii = 0; ii < L; ii++)
        {
            LTFAT_REAL inAbs = ltfat_abs(in[ii]);
            if (inAbs > normfac)
                normfac = inAbs;
        }
        break;

    }
    case LTFAT_NORMALIZE_NULL:
        normfac = 1.0;
        break;
    default:
        CHECKCANTHAPPEN("Unknown normalization flag");
    };

    normfac = (LTFAT_REAL)(1.0) / normfac;

    for (ltfat_int ii = 0; ii < L; ii++)
        out[ii] = normfac * in[ii];

error:
    return status;
}

/* LTFAT_API int */
/* LTFAT_NAME(postpad)(const LTFAT_TYPE* in, ltfat_int Ls, ltfat_int W, */
/*                     ltfat_int L, LTFAT_TYPE* out) */
/* { */
/*     int status = LTFATERR_SUCCESS; */
/*     CHECKNULL(in); CHECKNULL(out); */
/*     CHECK(LTFATERR_NOTPOSARG, Ls > 0, "Ls must be positive"); */
/*     CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive"); */
/*     CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive"); */
/*  */
/*     if (in == out) */
/*     { */
/*         LTFAT_TYPE* outTmp = ltfat_malloc(L * W * sizeof * out); */
/*         CHECKMEM(outTmp); */
/*     } */
/*     else */
/*     { */
/*         outTmp = out; */
/*     } */
/*  */
/*     ltfat_int Lcom = (Ls < L ? Ls : L); */
/*     ltfat_int Lrem = L - Lcom; */
/*  */
/*     for (ltfat_int w = 0; w < W; w++) */
/*     { */
/*         memcpy(outTmp + w * L, in + w * Ls, Lcom * sizeof * out); */
/*         memset(outTmp + w * L + Ls, 0, Lrem * sizeof * out); */
/*     } */
/*  */
/*     if (in == out) */
/*     { */
/*         ltfat_free(in); */
/*         out = outTmp; */
/*     } */
/*  */
/* error: */
/*     return status; */
/* } */



/* LTFAT_API LTFAT_REAL */
/* LTFAT_NAME(norm)(const LTFAT_TYPE in[], ltfat_int L, */
/*                  ltfat_normalize_t flag) */
/* { */
/*     double retNorm = 0.0; */
/*  */
/*     switch (flag) */
/*     { */
/*     case LTFAT_NORMALIZE_ENERGY: */
/*     { */
/*         for (ltfat_int ii = 0; ii < L; ii++) */
/*         { */
/* #ifdef LTFAT_COMPLEXTYPE */
/*             double inTmp = fabs(in[ii]); */
/*             retNorm += in[ii] * in[ii]; */
/* #else */
/*             retNorm += in[ii] * in[ii]; */
/* #endif */
/*         } */
/*     } */
/*     }; */
/*  */
/*     return (LTFAT_REAL) sqrt(retNorm); */
/* } */
