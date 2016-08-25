#ifndef _LTFAT_DST_H
#define _LTFAT_DST_H

typedef enum
{
    DSTI = FFTW_RODFT00, DSTIII = FFTW_RODFT01,
    DSTII = FFTW_RODFT10, DSTIV = FFTW_RODFT11
} dst_kind;

#endif

LTFAT_API LTFAT_FFTW(plan)
LTFAT_NAME(dst_init)( ltfat_int L, ltfat_int W, LTFAT_TYPE *cout,
                      const dst_kind kind);

LTFAT_API void
LTFAT_NAME(dst)(const LTFAT_TYPE *f, ltfat_int L, ltfat_int W,
                LTFAT_TYPE *cout, const dst_kind kind);

LTFAT_API void
LTFAT_NAME(dst_execute)(LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
                        ltfat_int L, ltfat_int W, LTFAT_TYPE *cout,
                        const dst_kind kind);
