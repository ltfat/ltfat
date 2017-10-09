#ifndef _ltfat_dgtrealmp_private_h
#define _ltfat_dgtrealmp_private_h
#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct ltfat_dgtrealmp_params
{
    ltfat_dgtrealmp_hint hint;
    ltfat_dgtrealmp_alg alg;
    long double errtoldb;
    long double errtoladj;
    double kernrelthr;
    size_t maxit;
    size_t maxatoms;
    size_t iterstep;
    int verbose;
    int initwasrun;
    int treelevels;
    ltfat_phaseconvention ptype;
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
    ltfat_int m;
    ltfat_int n;
    ltfat_int w;
    ltfat_int n2;
} kpoint;
#define PTOI(k) k.w][k.m + p->M2[k.w] * k.n
#define kpoint_init(m,n,w) LTFAT_STRUCTINIT(kpoint,m,n,w,n)
#define kpoint_init2(m,n,n2,w) LTFAT_STRUCTINIT(kpoint,m,n,w,n2)


typedef struct
{
    ksize           size;
    kanchor          mid;
    ltfat_int        kNo;
    LTFAT_COMPLEX** mods;
    LTFAT_COMPLEX** kval;
    double Mrat;
    double arat;
    ltfat_int Mstep;
    ltfat_int astep;
} LTFAT_NAME(kerns);


typedef struct
{
    LTFAT_REAL**           s;
    LTFAT_COMPLEX**        c;
    LTFAT_REAL**           maxcols;
    ltfat_int**            maxcolspos;
    LTFAT_NAME(maxtree)**  tmaxtree;
    LTFAT_NAME(maxtree)*** fmaxtree;
    LTFAT_REAL**           suppind;
    long double            err;
    long double            fnorm2;
    size_t                 currit;
    size_t                 curratoms;
    ltfat_int              P;
    ltfat_int*             N;
    // LocOMP related
    LTFAT_COMPLEX*         gramBuf;
    LTFAT_COMPLEX*         cvalBuf;
    LTFAT_COMPLEX*         cvalinvBuf;
    kpoint*                cvalBufPos;
    LTFAT_NAME_COMPLEX(hermsystemsolver_plan)* hplan;
    // CyclicMP related
    kpoint*                pointBuf;
    size_t                 pointBufNo;
} LTFAT_NAME(dgtrealmpiter_state);

// typedef struct
// {
//
// } LTFAT_NAME(dgtrealmpiter_atompos);


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
                                   int do_allmods, ltfat_phaseconvention ptype,
                                   LTFAT_NAME(kerns)** pout);

int
LTFAT_NAME(dgtrealmp_kernel_done)(LTFAT_NAME(kerns)** k);

int
LTFAT_NAME(dgtrealmp_kernel_modfi)(LTFAT_COMPLEX* kfirst, ksize size,
                                   kanchor mid, ltfat_int n, ltfat_int a, ltfat_int M,
                                   LTFAT_COMPLEX* kmod);
int
LTFAT_NAME(dgtrealmp_kernel_modti)(LTFAT_COMPLEX* kfirst, ksize size,
                                   kanchor mid, ltfat_int m, ltfat_int a, ltfat_int M,
                                   LTFAT_COMPLEX* kmod);

int
LTFAT_NAME(dgtrealmp_kernel_findsmallsize)(
    const LTFAT_COMPLEX kernlarge[], ltfat_int M, ltfat_int N,
    LTFAT_REAL reltol, ksize* size, kanchor* anchor);

int
LTFAT_NAME(dgtrealmp_essentialsupport)(
    const LTFAT_REAL g[], ltfat_int gl, LTFAT_REAL reltol,
    ltfat_int* lefttail, ltfat_int* righttail);

int
LTFAT_NAME(dgtrealmp_execute_kpos)(
    LTFAT_NAME(dgtrealmp_plan)* p, kpoint pos1, kpoint pos2,
    ltfat_int* m2, ltfat_int* n2, ltfat_int* Mstep, ltfat_int* astep,
    ksize* kdim2, kanchor* kmid2, kpoint* kstart2);

int
LTFAT_NAME(dgtrealmp_execute_indices)(
    LTFAT_NAME(dgtrealmp_plan)* p, kpoint origpos, kpoint* pos,
    ltfat_int* m2start, ltfat_int* n2start, ksize* kdim2, kanchor* kmid2,
    kpoint* kstart2);

LTFAT_COMPLEX*
LTFAT_NAME(dgtrealmp_execute_pickkernel)(
    LTFAT_NAME(kerns)* currkern, ltfat_int m, ltfat_int n,
    ltfat_phaseconvention pconv);

int
LTFAT_NAME(dgtrealmp_execute_findmaxatom)(
    LTFAT_NAME(dgtrealmp_plan)* p, kpoint* pos);
// ltfat_int* m, ltfat_int* n, ltfat_int* w);

int
LTFAT_NAME(dgtrealmp_execute_updateresiduum)(
    LTFAT_NAME(dgtrealmp_plan)* p, kpoint pos, LTFAT_COMPLEX cval,
    int do_substract);

LTFAT_REAL
LTFAT_NAME(dgtrealmp_execute_normfac)(
    LTFAT_NAME(dgtrealmp_plan)* p, kpoint pos, LTFAT_COMPLEX cval);

LTFAT_REAL
LTFAT_NAME(dgtrealmp_execute_mp)( LTFAT_NAME(dgtrealmp_plan)* p,
                                  kpoint pos, LTFAT_COMPLEX** cout);


LTFAT_REAL
LTFAT_NAME(dgtrealmp_execute_invmp)( LTFAT_NAME(dgtrealmp_plan)* p,
                                     kpoint pos, LTFAT_COMPLEX** cout);

// inline int
// LTFAT_NAME(dgtrealmp_execute_substractatom)();


#endif
