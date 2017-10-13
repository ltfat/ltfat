#ifndef _LTFAT_DGTREALMP_H
#define _LTFAT_DGTREALMP_H

enum ltfat_dgtrealmp_status
{
    LTFAT_DGTREALMP_STATUS_TOLREACHED   = 0,
    LTFAT_DGTREALMP_STATUS_MAXATOMS     = 1,
    LTFAT_DGTREALMP_STATUS_MAXITER      = 2,
    LTFAT_DGTREALMP_STATUS_STALLED      = 3,
    LTFAT_DGTREALMP_STATUS_LOCOMP_NOTHERM = 4,
    LTFAT_DGTREALMP_STATUS_LOCOMP_ORTHFAILED = 5,
    LTFAT_DGTREALMP_STATUS_CANCONTINUE  = 10
};


// typedef enum
// {
//     ltfat_dgtrealmp_auto         = 0,
//     ltfat_dgtrealmp_singlemod    = (1U << 0),
//     ltfat_dgtrealmp_allmods      = (1U << 1),
// } ltfat_dgtrealmp_hint;

typedef enum
{
    ltfat_dgtrealmp_alg_MP              = 0,
    ltfat_dgtrealmp_alg_LocOMP          = (1U << 0),
    ltfat_dgtrealmp_alg_LocCyclicMP     = (1U << 1),
} ltfat_dgtrealmp_alg;


typedef struct ltfat_dgtrealmp_params ltfat_dgtrealmp_params;

LTFAT_API ltfat_dgtrealmp_params*
ltfat_dgtrealmp_params_allocdef();

LTFAT_API int
ltfat_dgtrealmp_params_free(ltfat_dgtrealmp_params* params);

LTFAT_API int
ltfat_dgtrealmp_setpar_phaseconv(
    ltfat_dgtrealmp_params* params, ltfat_phaseconvention pconv);

// LTFAT_API int
// ltfat_dgtrealmp_setpar_hint(
//     ltfat_dgtrealmp_params* params, ltfat_dgtrealmp_hint hint);

LTFAT_API int
ltfat_dgtrealmp_setpar_alg(
    ltfat_dgtrealmp_params* params, ltfat_dgtrealmp_alg alg);

LTFAT_API int
ltfat_dgtrealmp_setpar_maxatoms(
    ltfat_dgtrealmp_params* params, size_t maxatoms);

LTFAT_API int
ltfat_dgtrealmp_setpar_maxit(
    ltfat_dgtrealmp_params* params, size_t maxit);

LTFAT_API int
ltfat_dgtrealmp_setpar_cycles(
    ltfat_dgtrealmp_params* params, size_t cycles);

LTFAT_API int
ltfat_dgtrealmp_setpar_kernrelthr(
    ltfat_dgtrealmp_params* p, double thr);

LTFAT_API int
ltfat_dgtrealmp_setpar_iterstep(
    ltfat_dgtrealmp_params* p, size_t iterstep);

LTFAT_API int
ltfat_dgtrealmp_setpar_errtoldb(
    ltfat_dgtrealmp_params* p, double errtoldb);

int
ltfat_dgtrealmp_params_defaults(ltfat_dgtrealmp_params* params);

// int
// ltfat_dgtrealmp_hint_isvalid(ltfat_dgtrealmp_hint in);

int
ltfat_dgtrealmp_alg_isvalid(ltfat_dgtrealmp_alg in);

#endif

typedef struct LTFAT_NAME(dgtrealmp_state) LTFAT_NAME(dgtrealmp_state);

LTFAT_API LTFAT_NAME(dgtreal_plan)**
LTFAT_NAME(dgtrealmp_getdgtrealplan)(LTFAT_NAME(dgtrealmp_state)* p);

LTFAT_API ltfat_dgtrealmp_params*
LTFAT_NAME(dgtrealmp_getparams)(const LTFAT_NAME(dgtrealmp_state)* p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_getresidualcoef_compact)(
    LTFAT_NAME(dgtrealmp_state)* p, LTFAT_COMPLEX* c);

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_iterstep)(
    LTFAT_NAME(dgtrealmp_state)* p, size_t iterstep);

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_maxatoms)(
    LTFAT_NAME(dgtrealmp_state)* p, size_t maxatoms);

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_errtoldb)(
    LTFAT_NAME(dgtrealmp_state)* p, double errtoldb);

LTFAT_API int
LTFAT_NAME(dgtrealmp_get_errdb)(
    const LTFAT_NAME(dgtrealmp_state)* p, double* err);

LTFAT_API int
LTFAT_NAME(dgtrealmp_init_compact)(
    const LTFAT_REAL g[], ltfat_int gl[], ltfat_int L, ltfat_int P,
    ltfat_int a[], ltfat_int M[], ltfat_dgtrealmp_params* params,
    LTFAT_NAME(dgtrealmp_state)** pout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_init)(
    const LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P,
    ltfat_int a[], ltfat_int M[], ltfat_dgtrealmp_params* params,
    LTFAT_NAME(dgtrealmp_state)** p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_compact)(
    LTFAT_NAME(dgtrealmp_state)* p, const LTFAT_REAL* f,
    LTFAT_COMPLEX* cout, LTFAT_REAL* fout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute)(
    LTFAT_NAME(dgtrealmp_state)* p, const LTFAT_REAL* f,
    LTFAT_COMPLEX** cout, LTFAT_REAL* fout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_reset)(
    LTFAT_NAME(dgtrealmp_state)* p, const LTFAT_REAL* f);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_niters)(
    LTFAT_NAME(dgtrealmp_state)* p, ltfat_int itno, LTFAT_COMPLEX** cout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_niters_compact)(
    LTFAT_NAME(dgtrealmp_state)* p, ltfat_int itno, LTFAT_COMPLEX* cout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_revert)(
    LTFAT_NAME(dgtrealmp_state)* p, LTFAT_COMPLEX** cout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_revert_compact)(
    LTFAT_NAME(dgtrealmp_state)* p, LTFAT_COMPLEX* cout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_done)(LTFAT_NAME(dgtrealmp_state)** p);
