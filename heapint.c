#include "ltfat.h"
#include "ltfat_types.h"
#include "heapint.h"

#define NORTHFROMW(w,M,N) ((((w) + 1) % (M)) + (w) - (w) % (M))
#define SOUTHFROMW(w,M,N) (((w) - 1 + (M)) % (M) + (w) - (w) % (M))

#define EASTFROMW(w,M,N)  (((w) + (M)) % ((M) * (N)))
#define WESTFROMW(w,M,N)  (((w) - (M) + (M) * (N)) % ((M) * (N)))

/*
 * Heap is a dynamic array h of heapsize elements which are ordered according
 * to s such that s[h[0]] is the maximum.
 *
 * */
struct LTFAT_NAME(heap)
{
    ltfatInt *h;
    ltfatInt heapsize;
    ltfatInt totalheapsize;
    const LTFAT_REAL *s;
};

struct LTFAT_NAME(heapinttask)
{
    ltfatInt M;
    ltfatInt N;
    const LTFAT_REAL *tgrad;
    const LTFAT_REAL *fgrad;
    int *donemask;
};


LTFAT_EXTERN
void LTFAT_NAME(heap_insert)(struct LTFAT_NAME(heap) *h, const ltfatInt key)
{
    ltfatInt pos, pos2, swap;

    /* Grow heap if necessary */
    if (h->totalheapsize == h->heapsize)
    {
        (h->totalheapsize) *= 2;
        h->h = ltfat_realloc_and_copy(h->h,
                                      h->totalheapsize * sizeof * h->h / 2,
                                      h->totalheapsize * sizeof * h->h);
    }

    pos = h->heapsize;
    h->heapsize++;

    h->h[h->heapsize - 1] = key;

    while (pos > 0)
    {
        /* printf("pos %i\n",pos); */
        pos2 = (pos - pos % 2) / 2;
        if (h->s[h->h[pos2]] < h->s[h->h[pos]] )
        {
            swap = h->h[pos2];
            h->h[pos2] = h->h[pos];
            h->h[pos] = swap;
            pos = pos2;
        }
        else
        {
            break;
        }
    }
}

LTFAT_EXTERN
ltfatInt LTFAT_NAME(heap_delete)(struct LTFAT_NAME(heap) *h)
{

    ltfatInt pos, maxchildpos, swap, key;
    LTFAT_REAL maxchildkey;

    /* Extract first element */
    key = h->h[0];

    /* Put last element on first elements place, and make the heap smaller. */
    h->h[0] = h->h[h->heapsize - 1];
    h->heapsize--;

    /* Fix the just introduced problem. */
    pos = 0;

    /*  %%%%%%%%%%%%%%
     * %
     * %  Is maxchildpos 0 or 1 indexed!
     * %
     * %
     * %%%%%%%%%%%%%
     */

    while (2 * pos + 1 < h->heapsize)
    {
        if (2 * pos + 3 > h->heapsize)
        {
            maxchildkey = h->s[h->h[2 * pos + 1]];
            maxchildpos = 1;
        }
        else
        {
            if (h->s[h->h[2 * pos + 1]] >= h->s[h->h[2 * pos + 2]])
            {
                maxchildkey = h->s[h->h[2 * pos + 1]];
                maxchildpos = 1;
            }
            else
            {
                maxchildkey = h->s[h->h[2 * pos + 2]];
                maxchildpos = 2;
            }
        }

        if (maxchildkey > h->s[h->h[pos]])
        {
            swap = h->h[2 * pos + maxchildpos];
            h->h[2 * pos + maxchildpos] = h->h[pos];
            h->h[pos] = swap;
            pos = 2 * pos + maxchildpos;
        }
        else
        {
            break;
        }
    }

    return key;
}


void LTFAT_NAME(trapezheap)(struct LTFAT_NAME(heap) *h,
                            const struct LTFAT_NAME(heapinttask) *heaptask,
                            const ltfatInt w,
                            LTFAT_REAL* phase)
{
    const ltfatInt M = heaptask->M;
    const ltfatInt N = heaptask->N;
    const LTFAT_REAL* tgradw = heaptask->tgrad;
    const LTFAT_REAL* fgradw = heaptask->fgrad;
    int* donemask = heaptask->donemask;
    ltfatInt w_E, w_W, w_N, w_S;

    /* Try and put the four neighbours onto the heap.
     * Integration by trapezoidal rule */

    /* North */
    w_N = NORTHFROMW(w,M,N);
    if (!donemask[w_N])
    {
        phase[w_N] = phase[w] + (fgradw[w] + fgradw[w_N]) / 2;
        donemask[w_N] = 1;
        LTFAT_NAME(heap_insert)(h, w_N);
    }

    /* South */
    w_S = SOUTHFROMW(w,M,N);
    if (!donemask[w_S])
    {
        phase[w_S] = phase[w] - (fgradw[w] + fgradw[w_S]) / 2;
        donemask[w_S] = 2;
        LTFAT_NAME(heap_insert)(h, w_S);
    }

    /* East */
    w_E = EASTFROMW(w,M,N);
    if (!donemask[w_E])
    {
        phase[w_E] = phase[w] + (tgradw[w] + tgradw[w_E]) / 2;
        donemask[w_E] = 3;
        LTFAT_NAME(heap_insert)(h, w_E);
    }

    /* West */
    w_W = WESTFROMW(w,M,N);
    if (!donemask[w_W])
    {
        phase[w_W] = phase[w] - (tgradw[w] + tgradw[w_W]) / 2;
        donemask[w_W] = 4;
        LTFAT_NAME(heap_insert)(h, w_W);
    }
}

LTFAT_EXTERN
void LTFAT_NAME(heapint)(const LTFAT_REAL *s,
                         const LTFAT_REAL *tgrad,
                         const LTFAT_REAL *fgrad,
                         const ltfatInt a, const ltfatInt M,
                         const ltfatInt L, const ltfatInt W,
                         LTFAT_REAL tol, LTFAT_REAL *phase)
{

    /* Declarations */
    ltfatInt N, ii, Imax, domainloop;
    ltfatInt w;
    LTFAT_REAL maxs;
    int *donemask;
    struct LTFAT_NAME(heap) h;
    struct LTFAT_NAME(heapinttask) hit;

    N = L / a;

    /* Main body */

    h.totalheapsize  = (ltfatInt)(M * log((LTFAT_REAL)M));
    h.h              = ltfat_malloc(h.totalheapsize * sizeof * h.h);
    h.s              = s;
    h.heapsize       = 0;

    donemask = ltfat_malloc(M * N * W * sizeof * donemask);

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL *tgradw = ltfat_malloc(M * N * sizeof * tgradw);
    LTFAT_REAL *fgradw = ltfat_malloc(M * N * sizeof * fgradw);

    /* Set the phase to zero initially */
    memset(phase, 0, M * N * W * sizeof * phase);

    /* Rescale the derivatives such that they are in readians and the step is 1 */
    LTFAT_NAME(gradsamptorad)(tgrad,fgrad, a, M, L, tgradw, fgradw);

    /* We will start intergration from the biffest coefficient */
    LTFAT_NAME_REAL(findmaxinarray)(s,M*N,&maxs,&Imax);

    /* Mark all the small elements as done, they get a zero phase.
     * Code 5
     */
    for (ii = 0; ii < M * N * W; ii++)
    {
        if (s[ii] < tol * maxs)
            donemask[ii] = 5;
        else
            donemask[ii] = 0;
    }

    /* Create a struct holding inputs */
    hit = (struct LTFAT_NAME(heapinttask)) { M, N, tgradw, fgradw, donemask };

    /* Outer loop over islands of coefficients */
    domainloop = 1;
    while (domainloop)
    {
        /* Empty the heap */
        h.heapsize = 0;

        /* Put maximal element onto the heap and mark that it is done. It
         * will get zero phase
         */
        LTFAT_NAME(heap_insert)(&h, Imax);
        donemask[Imax] = 6;

        /* Inner loop processing all connected coefficients */
        while (h.heapsize > 0)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(&h);
            /* Spread the current phase value to 4 direct neighbors */
            LTFAT_NAME(trapezheap)(&h, &hit, w, phase);
        }

        /* Find the new maximal element */
        /* Break from the outer loop when there is none */
        domainloop = LTFAT_NAME_REAL(findmaxinarraywrtmask)(s,donemask,M*N,&maxs,&Imax);
    }

    ltfat_free(tgradw);
    ltfat_free(fgradw);
    ltfat_free(donemask);
    ltfat_free(h.h);
}


    /* Rescale the partial derivatives to a torus of length L. We copy
       to work arrays because the input is const. */
    /* This is a 3 step process:
     *
     * 1) fgrad only: convert to frequency invariant
     * 2) from samples to radians
     * 3) Change step to 1
     *
     * */

void
LTFAT_NAME(gradsamptorad)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                          ltfatInt a, ltfatInt M, ltfatInt L,
                          LTFAT_REAL* tgradw, LTFAT_REAL* fgradw)
{
    ltfatInt N = L/a;
    ltfatInt b = L/M;
    LTFAT_REAL sampToRadConst = (LTFAT_REAL)( 2.0 * PI / L);

    for (ltfatInt jj = 0; jj < N; jj++)
    {
        for (ltfatInt ii = 0; ii < M; ii++)
        {
            tgradw[ii + jj * M] =    a * tgrad[ii + jj * M] * sampToRadConst;
            fgradw[ii + jj * M] =  - b * ( fgrad[ii + jj * M] + jj * a ) * sampToRadConst;
            /* The following converts phase derivatives so that the result is time-invariant phase */
            /* tgradw[ii + jj * M] =    a * (tgrad[ii + jj * M] + ii*b) * sampToRadConst; */
            /* fgradw[ii + jj * M] =  - b * ( fgrad[ii + jj * M] ) * sampToRadConst; */
        }
    }
}




#undef NORTHFROMW
#undef SOUTHFROMW
#undef WESTFROMW
#undef EASTFROMW

#undef NORTHEASTFROMW
#undef NORTHWESTFROMW
#undef SOUTHEASTFROMW
#undef SOUTHWESTFROMW
