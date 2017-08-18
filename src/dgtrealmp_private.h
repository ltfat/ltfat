#ifndef _ltfat_dgtrealmp_private_h
#define _ltfat_dgtrealmp_private_h
#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct ltfat_dgtrealmp_params
{
    ltfat_dgtrealmp_hint hint;
    ltfat_dgtrealmp_alg alg;
    LTFAT_REAL errtoldb;
    LTFAT_REAL errtoladj;
    LTFAT_REAL kernrelthr;
    size_t maxit;
    size_t maxatoms;
    size_t iterstep;
    int verbose;
    int initwasrun;
};

typedef struct
{
    ltfat_int height;
    ltfat_int width;
} ksize;

typedef struct
{
    ltfat_int hmid;
    ltfat_int wmid;
} kanchor;

typedef struct
{
    ksize        size;
    kanchor       mid;
    ltfat_int     kNo;
    LTFAT_COMPLEX** kval;
} LTFAT_NAME(kerns);

typedef struct
{
    LTFAT_REAL** s;
    LTFAT_COMPLEX** c;
    int** suppindCount;
    LTFAT_REAL** maxcols;
    ltfat_int** maxcolspos;
    LTFAT_REAL err;
    LTFAT_REAL fnorm2;
    size_t currit;
    size_t curratoms;
    ltfat_int         P;
} LTFAT_NAME(dgtrealmpiter_state);

struct LTFAT_NAME(dgtrealmp_plan)
{
    LTFAT_NAME(dgtrealmpiter_state)* iterstate;
    LTFAT_NAME(kerns)**             gramkerns; // PxP plans
    LTFAT_NAME(dgtreal_plan)**       dgtplans;  // P plans
    ltfat_int*        a;
    ltfat_int*        M;
    ltfat_int*       M2;
    ltfat_int*        N;
    ltfat_int         P;
    ltfat_int         L;
    LTFAT_REAL**    cout;
    ltfat_dgtrealmp_params* params;
};

int
LTFAT_NAME(dgtrealmpiter_init)(ltfat_int a[], ltfat_int M[], ltfat_int P,
                               ltfat_int L, LTFAT_NAME(dgtrealmpiter_state)** state);

int
LTFAT_NAME(dgtrealmpiter_done)(LTFAT_NAME(dgtrealmpiter_state)** state);

int
LTFAT_NAME(dgtrealmp_kernel_init)( const LTFAT_REAL* g[], ltfat_int gl[],
                                   ltfat_int a[], ltfat_int M[],
                                   ltfat_int L, LTFAT_REAL reltol,
                                   LTFAT_NAME(kerns)** pout);

int
LTFAT_NAME(dgtrealmp_kernel_done)(LTFAT_NAME(kerns)** k);

int
LTFAT_NAME(dgtrealmp_kernel_modfi)(LTFAT_COMPLEX* kfirst, ksize size,
                                   kanchor mid, ltfat_int n, ltfat_int a, ltfat_int M,
                                   LTFAT_COMPLEX* kmod);
int
LTFAT_NAME(dgtrealmp_kernel_findsmallsize)(const LTFAT_COMPLEX kernlarge[],
        ltfat_int M, ltfat_int N, LTFAT_REAL reltol,
        ksize* size, kanchor* anchor);

int
LTFAT_NAME(dgtrealmp_essentialsupport)(const LTFAT_REAL g[], ltfat_int gl,
                                       LTFAT_REAL reltol,
                                       ltfat_int* lefttail, ltfat_int* righttail);

#endif
