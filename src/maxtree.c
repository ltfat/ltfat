#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

int
LTFAT_NAME(maxtree_updatedirty)(LTFAT_NAME(maxtree)* p);

int
LTFAT_NAME(maxtree_updaterange)(
    LTFAT_NAME(maxtree)* p, ltfat_int start, ltfat_int stop);

struct LTFAT_NAME(maxtree)
{
    ltfat_int dirtystart;
    ltfat_int dirtyend;
    LTFAT_REAL*  pointedarray;
    LTFAT_REAL*  treeVals;
    LTFAT_REAL** treePtrs;
    ltfat_int*   treePos;
    ltfat_int**  treePosPtrs;
    ltfat_int    depth;
    ltfat_int    L;
    ltfat_int    Lstep;
    ltfat_int    nextL;
    ltfat_int*   levelL;
    ltfat_int    W;
};

LTFAT_API int
LTFAT_NAME(maxtree_init)(
    ltfat_int L, ltfat_int Lstep, ltfat_int depth,
    LTFAT_NAME(maxtree)** pout)
{
    LTFAT_NAME(maxtree)* p = NULL;
    ltfat_int nextL, granL, cumL;
    int status = LTFATERR_SUCCESS;

    CHECK(LTFATERR_NOTPOSARG, L > 0,
          "L must be positive (passed %td)" , L);
    CHECK(LTFATERR_BADARG, depth >= 0,
          "depth must be zero or greater (passed %td)", depth);

    granL = 1 << (depth);

    while ( granL > 2 * L )
        granL = 1 << (--depth);

    nextL = granL * ((L + granL - 1 ) / granL);

    CHECKMEM( p = LTFAT_NEW( LTFAT_NAME(maxtree)) );
    CHECKMEM( p->levelL = LTFAT_NEWARRAY(ltfat_int, depth + 1) );
    CHECKMEM( p->treePtrs = LTFAT_NEWARRAY(LTFAT_REAL*, depth + 1) );

    for (ltfat_int d = 0; d < depth + 1; d++)
    {
        ltfat_int dpow2 = 1 << d;
        p->levelL[depth - d] = (L + dpow2 - 1) / dpow2;
    }

    if (depth > 0)
    {

        CHECKMEM( p->treeVals = LTFAT_NAME_REAL(calloc)( nextL));
        CHECKMEM( p->treePos = LTFAT_NEWARRAY(ltfat_int, nextL ) );
        CHECKMEM( p->treePosPtrs = LTFAT_NEWARRAY(ltfat_int*, depth ) );

        cumL = 0;
        for (ltfat_int d = 0; d < depth; d++)
        {
            p->treePosPtrs[d] = p->treePos + cumL;
            p->treePtrs[d] = p->treeVals + cumL;
            cumL += p->levelL[d];
        }
    }

    p->depth = depth; p->L = L; p->nextL = nextL; p->Lstep = Lstep;

    p->dirtystart = p->Lstep;
    p->dirtyend   = 0;

    *pout = p;
    return LTFATERR_SUCCESS;
error:
    if (p) LTFAT_NAME(maxtree_done)(&p);
    return status;
}

LTFAT_API int
LTFAT_NAME(maxtree_done)(LTFAT_NAME(maxtree)** p)
{
    LTFAT_NAME(maxtree)* pp = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    ltfat_safefree(pp->treeVals);
    ltfat_safefree(pp->treePos);
    ltfat_safefree(pp->levelL);
    ltfat_safefree(pp->treePtrs);
    ltfat_safefree(pp->treePosPtrs);

    ltfat_free(pp);
    *p = NULL;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(maxtree_initwitharray)(
    ltfat_int L, ltfat_int depth, const LTFAT_REAL inarray[],
    LTFAT_NAME(maxtree)** pout)
{
    LTFAT_NAME(maxtree)* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS( LTFAT_NAME(maxtree_init)( L, L, depth, &p),
                 "init failed");

    CHECKSTATUS( LTFAT_NAME(maxtree_reset)( p, inarray),
                 "reset failed");

    *pout = p;
    return LTFATERR_SUCCESS;
error:
    if (p) LTFAT_NAME(maxtree_done)(&p);
    return status;
}


LTFAT_API int
LTFAT_NAME(maxtree_reset)(
    LTFAT_NAME(maxtree)* p, const LTFAT_REAL inarray[])
{
    p->treePtrs[p->depth] = (LTFAT_REAL*) inarray;

    return LTFAT_NAME(maxtree_updaterange)(p, 0, p->L);
}


LTFAT_API int
LTFAT_NAME(maxtree_setdirty)(LTFAT_NAME(maxtree)* p, ltfat_int start,
                             ltfat_int end)
{
    if (start < p->dirtystart) p->dirtystart = start;
    if (end   > p->dirtyend)  p->dirtyend  = end;
    return 0;
}

LTFAT_API int
LTFAT_NAME(maxtree_getdirty)(LTFAT_NAME(maxtree)* p, ltfat_int* start,
                             ltfat_int* end)
{
    *start = p->dirtystart;
    *end   = p->dirtyend;
    return 0;
}

int
LTFAT_NAME(maxtree_updatedirty)(LTFAT_NAME(maxtree)* p)
{
    if ( p->dirtyend <= p->dirtystart )
        return 1;

    int ret = LTFAT_NAME(maxtree_updaterange)(p, p->dirtystart, p->dirtyend);
    p->dirtystart = p->Lstep;
    p->dirtyend   = 0;
    return ret;
}

int
LTFAT_NAME(maxtree_updaterange)(LTFAT_NAME(maxtree)* p, ltfat_int start,
                                ltfat_int end)
{
    if (p->depth == 0) return 0;

    if (end > p->Lstep)
    {
        ltfat_int over = end - p->Lstep;
        LTFAT_NAME(maxtree_updaterange)( p, 0, over);
    }

    if (end > p->L) end = p->L;
    if (start >= end) return 0;

    ltfat_int parity = 0;
    parity =  end == p->L ? end % 2 : 0;
    start = start - start % 2;
    end   = end   + end % 2;
    start /= 2; end /= 2;

    LTFAT_REAL* treeVal = p->treePtrs[p->depth];
    LTFAT_REAL* treeValnext = p->treePtrs[p->depth - 1];
    ltfat_int* treePosnext = p->treePosPtrs[p->depth - 1];

    for (ltfat_int l = start; l < end - parity; l++)
    {
        if ( treeVal[2 * l] > treeVal[2 * l + 1])
        {
            treeValnext[l] = treeVal[2 * l];
            treePosnext[l] = 2 * l;
        }
        else
        {
            treeValnext[l] = treeVal[2 * l + 1];
            treePosnext[l] = 2 * l + 1;
        }
    }

    if ( parity )
    {
        treeValnext[end - 1] = treeVal[2 * (end - 1)];
        treePosnext[end - 1] = 2 * (end - 1);
    }

    for (ltfat_int d = p->depth - 1; d > 0; d--)
    {
        parity =  end >= p->levelL[d] ? end % 2 : 0;
        start = start - start % 2;
        end   = end   + end % 2;
        start /= 2; end /= 2;

        LTFAT_REAL* treeVal = p->treePtrs[d];
        LTFAT_REAL* treeValnext = p->treePtrs[d - 1];
        ltfat_int* treePos = p->treePosPtrs[d];
        ltfat_int* treePosnext = p->treePosPtrs[d - 1];

        for (ltfat_int l = start; l < end - parity; l++)
        {
            if ( treeVal[2 * l] > treeVal[2 * l + 1])
            {
                treeValnext[l] = treeVal[2 * l];
                treePosnext[l] = treePos[2 * l];
            }
            else
            {
                treeValnext[l] = treeVal[2 * l + 1];
                treePosnext[l] = treePos[2 * l + 1];
            }
        }

        if ( parity )
        {
            treeValnext[end - 1] = treeVal[2 * (end - 1)];
            treePosnext[end - 1] = treePos[2 * (end - 1)];
        }

    }

    return 0;
}

LTFAT_API int
LTFAT_NAME(maxtree_findmax)(LTFAT_NAME(maxtree)* p, LTFAT_REAL* max,
                            ltfat_int* maxPos)
{
    LTFAT_NAME(maxtree_updatedirty)(p);

    LTFAT_NAME_REAL(findmaxinarray)(p->treePtrs[0], p->levelL[0],
                                    max, maxPos);

    if (p->depth > 0)
        *maxPos = p->treePos[*maxPos];
    return 0;
}
