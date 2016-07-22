#ifndef _DCT_H
#define _DCT_H

typedef enum
{
    DCTI = FFTW_REDFT00, DCTIII = FFTW_REDFT01,
    DCTII = FFTW_REDFT10, DCTIV = FFTW_REDFT11
} dct_kind;

#endif

LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(dct_init)( const ltfatInt L, const ltfatInt W, LTFAT_TYPE *cout,
                      const dct_kind kind);


LTFAT_EXTERN void
LTFAT_NAME(dct)(const LTFAT_TYPE *f, const ltfatInt L, const ltfatInt W,
                LTFAT_TYPE *cout, const dct_kind kind);

LTFAT_EXTERN void
LTFAT_NAME(dct_execute)(const LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
                        const ltfatInt L, const ltfatInt W,
                        LTFAT_TYPE *cout, const dct_kind kind);
