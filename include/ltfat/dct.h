#ifndef _LTFAT_DCT_H
#define _LTFAT_DCT_H

typedef enum
{
    DCTI = FFTW_REDFT00, DCTIII = FFTW_REDFT01,
    DCTII = FFTW_REDFT10, DCTIV = FFTW_REDFT11
} dct_kind;

#endif

LTFAT_API LTFAT_FFTW(plan)
LTFAT_NAME(dct_init)( ltfat_int L, ltfat_int W, LTFAT_TYPE *cout,
                      const dct_kind kind);


LTFAT_API void
LTFAT_NAME(dct)(const LTFAT_TYPE *f, ltfat_int L, ltfat_int W,
                LTFAT_TYPE *cout, const dct_kind kind);

LTFAT_API void
LTFAT_NAME(dct_execute)(const LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
                        ltfat_int L, ltfat_int W,
                        LTFAT_TYPE *cout, const dct_kind kind);
