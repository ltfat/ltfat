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
    int treelevels;
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

struct LTFAT_NAME(maxtree)
{
    LTFAT_REAL* pointedarray;
    LTFAT_REAL* treeVals;
    LTFAT_REAL** treePtrs;
    ltfat_int*  treePos;
    ltfat_int** treePosPtrs;
    ltfat_int   depth;
    ltfat_int   L;
    ltfat_int   Lstep;
    ltfat_int   nextL;
    ltfat_int* levelL;
    ltfat_int   W;
};

typedef struct
{
    ksize        size;
    kanchor       mid;
    ltfat_int     kNo;
    LTFAT_COMPLEX** mods;
    LTFAT_COMPLEX** kval;
} LTFAT_NAME(kerns);

typedef struct
{
    LTFAT_REAL** s;
    LTFAT_COMPLEX** c;
    LTFAT_REAL** maxcols;
    ltfat_int**  maxcolspos;
    LTFAT_NAME(maxtree)** tmaxtree;
    LTFAT_NAME(maxtree)*** fmaxtree;
    int** suppindCount;
    LTFAT_REAL err;
    LTFAT_REAL fnorm2;
    size_t currit;
    size_t curratoms;
    ltfat_int         P;
    ltfat_int*        N;
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
    LTFAT_COMPLEX** couttmp;
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
                                   int do_allmods,
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
