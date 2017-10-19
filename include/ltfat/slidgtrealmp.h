#ifndef _LTFAT_SLIDGTREALMP_H
#define _LTFAT_SLIDGTREALMP_H


#endif

typedef struct LTFAT_NAME(slidgtrealmp_state) LTFAT_NAME(slidgtrealmp_state);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_init)(
    LTFAT_NAME(dgtrealmp_state)* p, int take_ownership,
    LTFAT_NAME(slidgtrealmp_state)** pout);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute_compact)(
    LTFAT_NAME(slidgtrealmp_state)* p, const LTFAT_REAL f[],
    LTFAT_COMPLEX cout[], LTFAT_REAL fout[]);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute)(
    LTFAT_NAME(slidgtrealmp_state)* p, const LTFAT_REAL f[],
    LTFAT_COMPLEX* cout[], LTFAT_REAL fout[]);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_reset)(
    LTFAT_NAME(slidgtrealmp_state)* p, const LTFAT_REAL f[]);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_done)(LTFAT_NAME(slidgtrealmp_state)** p);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute_niters)(
    LTFAT_NAME(slidgtrealmp_state)* p, ltfat_int itno, LTFAT_COMPLEX* cout[]);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute_niters_compact)(
    LTFAT_NAME(slidgtrealmp_state)* p, ltfat_int itno, LTFAT_COMPLEX cout[]);

