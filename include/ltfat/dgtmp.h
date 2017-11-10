#ifndef _LTFAT_DGTMP_H
#define _LTFAT_DGTMP_H

enum ltfat_dgtmp_status
{
    LTFAT_DGTMP_STATUS_TOLREACHED   = 0,
    LTFAT_DGTMP_STATUS_MAXATOMS     = 1,
    LTFAT_DGTMP_STATUS_MAXITER      = 2,
    LTFAT_DGTMP_STATUS_STALLED      = 3,
    LTFAT_DGTMP_STATUS_LOCOMP_NOTHERM = 4,
    LTFAT_DGTMP_STATUS_LOCOMP_ORTHFAILED = 5,
    LTFAT_DGTMP_STATUS_CANCONTINUE  = 10
};

typedef enum
{
    ltfat_dgtmp_alg_MP              = 0,
    ltfat_dgtmp_alg_LocOMP          = 1,
    ltfat_dgtmp_alg_LocCyclicMP     = 2,
} ltfat_dgtmp_alg;

typedef struct ltfat_dgtmp_params ltfat_dgtmp_params;

LTFAT_API ltfat_dgtmp_params*
ltfat_dgtmp_params_allocdef();

LTFAT_API int
ltfat_dgtmp_params_free(ltfat_dgtmp_params* params);

LTFAT_API int
ltfat_dgtmp_setpar_phaseconv(
    ltfat_dgtmp_params* params, ltfat_phaseconvention pconv);

// LTFAT_API int
// ltfat_dgtmp_setpar_hint(
//     ltfat_dgtmp_params* params, ltfat_dgtmp_hint hint);

LTFAT_API int
ltfat_dgtmp_setpar_alg(
    ltfat_dgtmp_params* params, ltfat_dgtmp_alg alg);

LTFAT_API int
ltfat_dgtmp_setpar_maxatoms(
    ltfat_dgtmp_params* params, size_t maxatoms);

LTFAT_API int
ltfat_dgtmp_setpar_maxit(
    ltfat_dgtmp_params* params, size_t maxit);

LTFAT_API int
ltfat_dgtmp_setpar_cycles(
    ltfat_dgtmp_params* params, size_t cycles);

LTFAT_API int
ltfat_dgtmp_setpar_kernrelthr(
    ltfat_dgtmp_params* p, double thr);

LTFAT_API int
ltfat_dgtmp_setpar_iterstep(
    ltfat_dgtmp_params* p, size_t iterstep);

LTFAT_API int
ltfat_dgtmp_setpar_errtoldb(
    ltfat_dgtmp_params* p, double errtoldb);

int
ltfat_dgtmp_params_defaults(ltfat_dgtmp_params* params);

// int
// ltfat_dgtmp_hint_isvalid(ltfat_dgtmp_hint in);

int
ltfat_dgtmp_alg_isvalid(ltfat_dgtmp_alg in);

#endif

typedef struct LTFAT_NAME(dgtmp_state) LTFAT_NAME(dgtmp_state);

LTFAT_API int
LTFAT_NAME(dgtmp_init)(
    const LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[],
    ltfat_int M[], ltfat_dgtmp_params* params, LTFAT_NAME(dgtmp_state)** pout);

LTFAT_API int
LTFAT_NAME(dgtmp_reset)(LTFAT_NAME(dgtmp_state)* p, const LTFAT_TYPE* f);

LTFAT_API int
LTFAT_NAME(dgtmp_execute_niters)(
    LTFAT_NAME(dgtmp_state)* p, ltfat_int itno, LTFAT_COMPLEX** cout);

LTFAT_API int
LTFAT_NAME(dgtmp_execute)(
    LTFAT_NAME(dgtmp_state)* p,
    const LTFAT_TYPE* f, LTFAT_COMPLEX** cout, LTFAT_COMPLEX* fout);

LTFAT_API int
LTFAT_NAME(dgtmp_done)( LTFAT_NAME(dgtmp_state)** p);

LTFAT_API int
LTFAT_NAME(dgtmp_init_compact)(
    const LTFAT_REAL g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[],
    ltfat_int M[], ltfat_dgtmp_params* params,
    LTFAT_NAME(dgtmp_state)** pout);

LTFAT_API int
LTFAT_NAME(dgtmp_execute_compact)(
    LTFAT_NAME(dgtmp_state)* p,
    const LTFAT_TYPE* f, LTFAT_COMPLEX* cout, LTFAT_COMPLEX* fout);


LTFAT_API int
LTFAT_NAME(dgtmp_execute_niters_compact)(
    LTFAT_NAME(dgtmp_state)* p, ltfat_int itno, LTFAT_COMPLEX* cout);

LTFAT_API int
LTFAT_NAME(dgtmp_revert)(
    LTFAT_NAME(dgtmp_state)* p, LTFAT_COMPLEX** cout);
