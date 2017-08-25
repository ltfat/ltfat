#ifndef _LTFAT_DGTREALMP_H
#define _LTFAT_DGTREALMP_H

enum ltfat_dgtrealmp_status
{
LTFAT_DGTREALMP_STATUS_TOLREACHED = 0,
LTFAT_DGTREALMP_STATUS_MAXATOMS   = 1,
LTFAT_DGTREALMP_STATUS_MAXITER    = 2,
LTFAT_DGTREALMP_STATUS_STALLED    = 3,
LTFAT_DGTREALMP_STATUS_WILLCONTINUE   = 10
};


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
LTFAT_NAME(dgtrealmp_set_errtoldb)(LTFAT_NAME(dgtrealmp_plan)* p, int errtoldb);


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


// maxtree

typedef struct LTFAT_NAME(maxtree) LTFAT_NAME(maxtree);

LTFAT_API int
LTFAT_NAME(maxtree_init)( ltfat_int L, ltfat_int Lstep, ltfat_int depth, LTFAT_NAME(maxtree)** p);

LTFAT_API int
LTFAT_NAME(maxtree_initwitharray)(
    ltfat_int L, ltfat_int depth,
    const LTFAT_REAL* inarray,
    LTFAT_NAME(maxtree)** p);

LTFAT_API int
LTFAT_NAME(maxtree_reset)(LTFAT_NAME(maxtree)* p, const LTFAT_REAL* inarray);

LTFAT_API int
LTFAT_NAME(maxtree_updaterange)(LTFAT_NAME(maxtree)* p, ltfat_int start, ltfat_int stop);

LTFAT_API int
LTFAT_NAME(maxtree_findmax)(LTFAT_NAME(maxtree)* p, LTFAT_REAL* max, ltfat_int* maxPos);

LTFAT_API int
LTFAT_NAME(maxtree_done)(LTFAT_NAME(maxtree)** p);


