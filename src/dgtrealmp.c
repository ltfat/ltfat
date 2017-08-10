#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

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
    LTFAT_REAL** kval;
} LTFAT_NAME(kerns);

typedef struct
{
    LTFAT_REAL** s;
    LTFAT_REAL** c;
    int* suppindCount;
    LTFAT_REAL** maxcols;
    ltfat_int** maxcolspos;
    LTFAT_REAL err;
    size_t currit;
    size_t curratoms;
} LTFAT_NAME(dgtrealmpiter_state);

struct LTFAT_NAME(dgtrealmp_plan)
{
    LTFAT_NAME(dgtrealmpiter_state) iterstate;
    LTFAT_NAME(kerns)               gramkerns;
    LTFAT_NAME(dgtreal_plan)**       dgtplans;
    ltfat_int*        a;
    ltfat_int*        M;
    ltfat_int*       M2;
    ltfat_int*        N;
    LTFAT_REAL      tol;
    LTFAT_REAL**    cout;
    DGTREALMP_HINT hint;
    LTFAT_REAL errtol;
    size_t maxit;
    size_t maxatoms;
};


LTFAT_API int
LTFAT_NAME(dgtrealmp_init)(
    LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[],
    ltfat_int M[], DGTREALMP_HINT hint, LTFAT_NAME(dgtrealmp_plan)** pout)
{
    int status = LTFATERR_FAILED;
    LTFAT_NAME(dgtrealmp_plan)* p = NULL;

    ltfat_int nextL = ltfat_dgtlengthmulti(L, P, a, M);

    CHECK(LTFATERR_BADTRALEN, L == nextL,
          "Next compatible transform length is %d (passed %d).", nextL, L);

    CHECKMEM( p = LTFAT_NEW(LTFAT_NAME(dgtrealmp_plan)) );
    CHECKMEM( p->dgtplans = LTFAT_NEWARRAY(LTFAT_NAME(dgtreal_plan)*, P) );

    for (ltfat_int p = 0; p < P; p++)
    {
        CHECKSTATUS(
            LTFAT_NAME(dgtreal_init_gen)(g[p], gl[p], g[p], gl[p], L, 1, a[p], M[p],
                                         NULL, NULL, NULL, &p->dgtplans[p] ),
            "dgtreal_init failed" );
    }

    return LTFATERR_SUCCESS;
error:
    if (p) LTFAT_NAME(dgtrealmp_done)(&p);
    return status;
}
