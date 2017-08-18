#ifndef _LTFAT_DGTREALMP_H
#define _LTFAT_DGTREALMP_H

typedef enum
{
    ltfat_dgtrealmp_singlemod    = 0,
    ltfat_dgtrealmp_allmods      = (1U << 0),
} ltfat_dgtrealmp_hint;

typedef enum
{
    ltfat_dgtrealmp_alg_mp       = 0,
    ltfat_dgtrealmp_alg_locomp   = (1U << 0),
    ltfat_dgtrealmp_alg_cyclicmp  = (1U << 1),
} ltfat_dgtrealmp_alg;


typedef struct ltfat_dgtrealmp_params ltfat_dgtrealmp_params;

LTFAT_API ltfat_dgtrealmp_params*
ltfat_dgtrealmp_params_allocdef();

LTFAT_API int
ltfat_dgtrealmp_params_free(ltfat_dgtrealmp_params* params);

LTFAT_API int
ltfat_dgtrealmp_params_free(ltfat_dgtrealmp_params* params);



int
ltfat_dgtrealmp_params_defaults(ltfat_dgtrealmp_params* params);

#endif

typedef struct LTFAT_NAME(dgtrealmp_plan) LTFAT_NAME(dgtrealmp_plan);

LTFAT_API LTFAT_NAME(dgtreal_plan)**
LTFAT_NAME(dgtrealmp_getdgtrealplan)(LTFAT_NAME(dgtrealmp_plan)* p);

LTFAT_API ltfat_dgtrealmp_params*
LTFAT_NAME(dgtrealmp_getparams)(LTFAT_NAME(dgtrealmp_plan)* p);

LTFAT_API LTFAT_COMPLEX**
LTFAT_NAME(dgtrealmp_getresidualcoef)(LTFAT_NAME(dgtrealmp_plan)* p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_iterstep)(LTFAT_NAME(dgtrealmp_plan)* p, size_t iterstep);

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_maxatoms)(LTFAT_NAME(dgtrealmp_plan)* p, size_t maxatoms);


LTFAT_API int
LTFAT_NAME(dgtrealmp_init)(
    const LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[], ltfat_int M[],
    ltfat_dgtrealmp_params* params, LTFAT_NAME(dgtrealmp_plan)** p);


LTFAT_API int
LTFAT_NAME(dgtrealmp_execute)(LTFAT_NAME(dgtrealmp_plan)* p,
                              const LTFAT_REAL* f,
                              LTFAT_COMPLEX** cout,
                              LTFAT_REAL* fout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_reset)(LTFAT_NAME(dgtrealmp_plan)* p, const LTFAT_REAL* f);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_niters)(LTFAT_NAME(dgtrealmp_plan)* p, ltfat_int itno,
                                     LTFAT_COMPLEX** cout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_done)(LTFAT_NAME(dgtrealmp_plan)** p);
