#ifndef _LTFAT_DGTREALMP_H
#define _LTFAT_DGTREALMP_H

enum ltfat_dgtmp_status
{
    LTFAT_DGTREALMP_STATUS_TOLREACHED   = 0,
    LTFAT_DGTREALMP_STATUS_MAXATOMS     = 1,
    LTFAT_DGTREALMP_STATUS_MAXITER      = 2,
    LTFAT_DGTREALMP_STATUS_STALLED      = 3,
    LTFAT_DGTREALMP_STATUS_EMPTY        = 4,
    LTFAT_DGTREALMP_STATUS_LOCOMP_NOTHERM = 5,
    LTFAT_DGTREALMP_STATUS_LOCOMP_ORTHFAILED = 6,
    LTFAT_DGTREALMP_STATUS_CANCONTINUE  = 10
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
ltfat_dgtmp_setpar_kernrelthr(
    ltfat_dgtmp_params* p, double thr);

LTFAT_API int
ltfat_dgtmp_setpar_iterstep(
    ltfat_dgtmp_params* p, size_t iterstep);

LTFAT_API int
ltfat_dgtmp_setpar_errtoldb(
    ltfat_dgtmp_params* p, double errtoldb);

LTFAT_API int
ltfat_dgtmp_setpar_snrdb(
    ltfat_dgtmp_params* params, double snrdb);

LTFAT_API int
ltfat_dgtmp_setpar_pedanticsearch(
    ltfat_dgtmp_params* params, int do_pedanticsearch);

LTFAT_API int
ltfat_dgtmp_setpar_cycles(
        ltfat_dgtmp_params* params, size_t cycles);

// LTFAT_API int
// ltfat_dgtmp_setpar_checkerreverynit(
//     ltfat_dgtmp_params* p, ltfat_int itstep, double errtoldb);

int
ltfat_dgtmp_params_defaults(ltfat_dgtmp_params* params);

// int
// ltfat_dgtmp_hint_isvalid(ltfat_dgtmp_hint in);

int
ltfat_dgtmp_alg_isvalid(ltfat_dgtmp_alg in);

#endif

typedef struct LTFAT_NAME(dgtrealmp_state) LTFAT_NAME(dgtrealmp_state);
typedef struct LTFAT_NAME(dgtrealmp_parbuf) LTFAT_NAME(dgtrealmp_parbuf);

LTFAT_API LTFAT_NAME(dgtreal_plan)**
LTFAT_NAME(dgtrealmp_getdgtrealplan)(LTFAT_NAME(dgtrealmp_state)* p);

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
LTFAT_NAME(dgtrealmp_get_numatoms)(
    const LTFAT_NAME(dgtrealmp_state)* p, size_t* atoms);

LTFAT_API int
LTFAT_NAME(dgtrealmp_get_numiters)(
    const LTFAT_NAME(dgtrealmp_state)* p, size_t* iters);

LTFAT_API int
LTFAT_NAME(dgtrealmp_init)(
    LTFAT_NAME(dgtrealmp_parbuf)* pb, ltfat_int L, LTFAT_NAME(dgtrealmp_state)** p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_init_gen_compact)(
    const LTFAT_REAL g[], ltfat_int gl[], ltfat_int L, ltfat_int P,
    ltfat_int a[], ltfat_int M[], ltfat_dgtmp_params* params,
    LTFAT_NAME(dgtrealmp_state)** pout);

LTFAT_API int
LTFAT_NAME(dgtrealmp_init_gen)(
    const LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P,
    ltfat_int a[], ltfat_int M[], ltfat_dgtmp_params* params,
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

/***********************************************************************/


LTFAT_API ltfat_int
LTFAT_NAME(dgtrealmp_parbuf_nextcompatlen)(
    LTFAT_NAME(dgtrealmp_parbuf)* p, ltfat_int L);

LTFAT_API ltfat_int
LTFAT_NAME(dgtrealmp_parbuf_nextcoefsize)(
    LTFAT_NAME(dgtrealmp_parbuf) * p, ltfat_int L, ltfat_int dictid);

LTFAT_API ltfat_int
LTFAT_NAME(dgtrealmp_getparbuf_dictno)( LTFAT_NAME(dgtrealmp_parbuf) * p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_init)( LTFAT_NAME(dgtrealmp_parbuf)** p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_done)( LTFAT_NAME(dgtrealmp_parbuf)** p);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(
        LTFAT_NAME(dgtrealmp_parbuf)* parbuf,
        LTFAT_FIRWIN win, ltfat_int gl, ltfat_int a, ltfat_int M);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_add_genwin)(
        LTFAT_NAME(dgtrealmp_parbuf)* parbuf,
        const LTFAT_REAL g[], ltfat_int gl, ltfat_int a, ltfat_int M);

// LTFAT_API int
// LTFAT_NAME(dgtrealmp_parbuf_add_gausswin)(
//         LTFAT_NAME(dgtrealmp_parbuf)* parbuf,
//         LTFAT_REAL tfr, ltfat_int a, ltfat_int M);
//
// LTFAT_API int
// LTFAT_NAME(dgtrealmp_parbuf_add_hermwin)(
//         LTFAT_NAME(dgtrealmp_parbuf)* parbuf,
//         ltfat_int order, LTFAT_REAL tfr, ltfat_int a, ltfat_int M);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_phaseconv)(
        LTFAT_NAME(dgtrealmp_parbuf)* parbuf, ltfat_phaseconvention pconv);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_maxatoms)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, size_t maxatoms);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_maxit)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, size_t maxit);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_errtoldb)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, double errtoldb);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_snrdb)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, double snrdb);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_kernrelthr)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, double thr);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_iterstep)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, size_t iterstep);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_mod_makelasttight)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf);

/* TODO:
LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_mod_chirpmod)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, LTFAT_REAL chirprate);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_mod_shiftby)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, LTFAT_REAL frac);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_mod_truncat)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, LTFAT_REAL relthr);

LTFAT_API int
LTFAT_NAME(dgtrealmp_parbuf_add_longwin)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf,
    LTFAT_LONGWIN win, LTFAT_REAL damp, ltfat_int a, ltfat_int M);

// Callback allows changing the values of energies of inner products which
// are used in the search for the maximally corellated atom with the residuum.
typedef int LTFAT_NAME(selectionmodcallback)(void* userdata,
        const LTFAT_COMPLEX* in, ltfat_int, ltfat_int, ltfat_int, LTFAT_REAL* out);

LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_selectiomodcallback)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, LTFAT_NAME(selectionmodcallback)* callback);
*/
LTFAT_API int
LTFAT_NAME(dgtrealmp_setparbuf_alg)(
    LTFAT_NAME(dgtrealmp_parbuf)* parbuf, ltfat_dgtmp_alg alg);


