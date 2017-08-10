#ifndef _LTFAT_DGTREALMP_H
#define _LTFAT_DGTREALMP_H

typedef enum
{
    DGTREALMP_SINGLEMOD = 0,
    DGTREALMP_ALLMODS = (1U << 0)
} DGTREALMP_HINT;

#endif

typedef struct LTFAT_NAME(dgtrealmp_plan) LTFAT_NAME(dgtrealmp_plan);

LTFAT_API int
LTFAT_NAME(dgtrealmp_init)(
    LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[], ltfat_int M[],
    DGTREALMP_HINT hint, LTFAT_NAME(dgtrealmp_plan)** p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_reset)();

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute)(LTFAT_NAME(dgtrealmp_plan)* p,
                              const LTFAT_REAL* f, LTFAT_COMPLEX** cout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_sig)(LTFAT_NAME(dgtrealmp_plan)* p,
                                  const LTFAT_REAL* f, LTFAT_REAL* fout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_done)(LTFAT_NAME(dgtrealmp_plan)** p);


LTFAT_API int
LTFAT_NAME(dgtrealmpiter_init)();

LTFAT_API int
LTFAT_NAME(dgtrealmpiter_execute)();

LTFAT_API int
LTFAT_NAME(dgtrealmpiter_done)();
