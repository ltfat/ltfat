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


int
LTFAT_NAME(dgtrealmp_computekernel)(
    LTFAT_REAL g[],
    ltfat_int a[],
    ltfat_int M[],
    ltfat_int L,
    LTFAT_NAME(kerns)** k)
{
    double arat
    ltfat_int modNo, amin, Mmax, lefttail0, righttail0, lefttail1, righttail1;
    LTFAT_REAL* g0tmp, * g1tmp;
    int status = LTFATERR_FAILED;

    amin = a[0] > a[1] ? a[0] : a[1];
    Mmax = M[0] > M[1] ? M[0] : M[1];

    LTFAT_NAME(dgtrealmp_essentialsupport)(g[0], gl[0], 1e-6, &lefttail0,
                                           &righttail0);
    LTFAT_NAME(dgtrealmp_essentialsupport)(g[1], gl[1], 1e-6, &lefttail1,
                                           &righttail1);

    ltfat_int Lshort =
        (lefttail0 > righttail0 ? 2 * lefttail0 : 2 * righttail0) +
        (lefttail1 > righttail1 ? 2 * lefttail1 : 2 * righttail1);

    Lshort = ltfat_dgtlength(Lshort > L ? L : Lshort, amin, Mmax);
    Nshort = Lshort / amin;

    CHECKMEM(g0tmp = LTFAT_NAME_REAL(malloc)(Lshort));
    CHECKMEM(g1tmp = LTFAT_NAME_REAL(malloc)(Lshort));
    CHECKMEM(kernlarge = LTFAT_NAME_COMPLEX(malloc)(Mmax * Nshort));

    LTFAT_NAME(middlepad)(g[0], gl[0], Lshort, g0tmp);
    LTFAT_NAME(middlepad)(g[1], gl[1], Lshort, g1tmp);

    LTFAT_NAME_REAL(dgt_long)(g0tmp, g1tmp, Lshort, 1, amin, Mmax, LTFAT_FREQINV,
                              kernlarge);

    if (a[0] > a[1])
    {
        modNo = ltfat_lcm(amin, Mmax) / a[0];
        arat = a[0] / a[1];
    }
    else
    {
        modNo = ltfat_lcm(amin, Mmax) / amin;
        arat = 1;
    }





    return LTFATERR_SUCCESS;
error:
    return status;
}

int
LTFAT_NAME(dgtrealmp_essentialsupport)(const LTFAT_REAL g[], ltfat_int gl,
                                       LTFAT_REAL reltol,
                                       ltfat_int* lefttail, ltfat_int* righttail)
{
    ltfat_int gl2 = gl / 2 + 1;
    LTFAT_REAL gmax, gthr;
    ltfat_int gmaxPos, lastPos;

    LTFAT_NAME_REAL(findmaxinarray)(g, gl, &gmax, &gmaxPos);
    gthr = gmax * reltol;

    *righttail = 0;
    *lefttail  = 0;

    for (ltfat_int l = 0; l < gl2; l++)
        if (g[l] > gthr)
            *righttail = l;

    for (ltfat_int l = gl - 1; l > gl2; l--)
        if (g[l] > gthr)
            *lefttail = gl - l;

    return 0
}


int
LTFAT_NAME(dgtrealmp_findsmallkernelsize)(const LTFAT_REAL kernlarge[], ltfat_int M, ltfat_int N,
                                          LTFAT_REAL reltol,
                                          ksize* size, kanchor* anchor)
{
}


