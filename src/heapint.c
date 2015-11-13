#include "ltfat.h"
#include "ltfat_types.h"
#include "heapint.h"

#define NORTHFROMW(w,M,N) ((((w) + 1) % (M)) + (w) - (w) % (M))
#define SOUTHFROMW(w,M,N) (((w) - 1 + (M)) % (M) + (w) - (w) % (M))

#define EASTFROMW(w,M,N)  (((w) + (M)) % ((M) * (N)))
#define WESTFROMW(w,M,N)  (((w) - (M) + (M) * (N)) % ((M) * (N)))




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
    w_N = NORTHFROMW(w, M, N);
    if (!donemask[w_N])
    {
        phase[w_N] = phase[w] + (fgradw[w] + fgradw[w_N]) / 2;
        donemask[w_N] = 1;
        LTFAT_NAME(heap_insert)(h, w_N);
    }

    /* South */
    w_S = SOUTHFROMW(w, M, N);
    if (!donemask[w_S])
    {
        phase[w_S] = phase[w] - (fgradw[w] + fgradw[w_S]) / 2;
        donemask[w_S] = 2;
        LTFAT_NAME(heap_insert)(h, w_S);
    }

    /* East */
    w_E = EASTFROMW(w, M, N);
    if (!donemask[w_E])
    {
        phase[w_E] = phase[w] + (tgradw[w] + tgradw[w_E]) / 2;
        donemask[w_E] = 3;
        LTFAT_NAME(heap_insert)(h, w_E);
    }

    /* West */
    w_W = WESTFROMW(w, M, N);
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
    ltfatInt N, ii, Imax;
    ltfatInt w;
    LTFAT_REAL maxs;
    int *donemask;
    struct LTFAT_NAME(heap) h;
    struct LTFAT_NAME(heapinttask) hit;

    N = L / a;

    /* Main body */

    h.totalheapsize  = (ltfatInt)(M * log((LTFAT_REAL)M));
    h.h              = ltfat_malloc(h.totalheapsize * sizeof * h.h);
    h.s              = (LTFAT_REAL*)s;
    h.heapsize       = 0;

    donemask = ltfat_malloc(M * N * W * sizeof * donemask);

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL *tgradw = ltfat_malloc(M * N * sizeof * tgradw);
    LTFAT_REAL *fgradw = ltfat_malloc(M * N * sizeof * fgradw);

    /* Set the phase to zero initially */
    memset(phase, 0, M * N * W * sizeof * phase);

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptorad)(tgrad, fgrad, a, M, L, tgradw, fgradw);

    /* We will start intergration from the biffest coefficient */
    LTFAT_NAME_REAL(findmaxinarray)(s, M * N, &maxs, &Imax);

    /* Mark all the small elements as done, they get zero phase.
     * (But should get random phase instead)
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
    hit = (struct LTFAT_NAME(heapinttask))
    {
        M, N, tgradw, fgradw, donemask
    };

    /* Outer loop over islands of coefficients */
    do
    {
        /* Empty the heap */
        h.heapsize = 0;

        /* Put maximal element onto the heap and mark that it is done.
         * It gets a zero phase.  */
        LTFAT_NAME(heap_insert)(&h, Imax);
        donemask[Imax] = 6;

        /* Inner loop processing all connected coefficients */
        while (h.heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(&h);
            /* Spread the current phase value to 4 direct neighbors */
            LTFAT_NAME(trapezheap)(&h, &hit, w, phase);
        }

        /* Find the new maximal element */
        /* Break from the outer loop when there is none */
    }
    while (LTFAT_NAME_REAL(findmaxinarraywrtmask)(s, donemask, M * N, &maxs, &Imax));

    LTFAT_SAFEFREEALL(tgradw, fgradw, donemask, h.h);
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
    ltfatInt N = L / a;
    ltfatInt b = L / M;
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

LTFAT_EXTERN
void LTFAT_NAME(maskedheapint)(const LTFAT_COMPLEX  *c,
                               const LTFAT_REAL *tgrad,
                               const LTFAT_REAL *fgrad,
                               const int* mask,
                               const ltfatInt a, const ltfatInt M,
                               const ltfatInt L, const ltfatInt W,
                               LTFAT_REAL tol, LTFAT_REAL *phase)
{
    /* Declarations */
    ltfatInt N, ii, Imax;
    ltfatInt w;
    LTFAT_REAL maxs;
    int *donemask;
    struct LTFAT_NAME(heap) h;
    struct LTFAT_NAME(heapinttask) hit;

    N = L / a;

    /* Main body */
    h.s              = ltfat_malloc(M * N * sizeof * h.s);
    h.totalheapsize  = (ltfatInt)(M * log((LTFAT_REAL)M));
    h.h              = ltfat_malloc(h.totalheapsize * sizeof * h.h);
    h.heapsize       = 0;

    for (ltfatInt w = 0; w < M * N * W; w++)
    {
        h.s[w] = LTFAT_COMPLEXH(cabs)(c[w]);
    }

    donemask = ltfat_malloc(M * N * W * sizeof * donemask);

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL *tgradw = ltfat_malloc(M * N * sizeof * tgradw);
    LTFAT_REAL *fgradw = ltfat_malloc(M * N * sizeof * fgradw);

    /* Set the phase to zero initially */
    // memset(phase, 0, M * N * W * sizeof * phase);

    /* Copy known phase */
    /* Code 12 */
    for (ltfatInt w = 0; w < M * N * W; w++)
    {
        if (mask[w])
        {
            phase[w]    = LTFAT_COMPLEXH(carg)(c[w]);
            donemask[w] = 12; /* Code of known phase */
        }
        else
        {
            phase[w]    = 0;
            donemask[w] = 0; /* Mask of unknown phase */
        }
    }

    LTFAT_NAME_REAL(findmaxinarray)(h.s, M * N, &maxs, &Imax);
    /* Mark all the small elements as done, they get zero phase.
     * Code 5
     */
    for (ii = 0; ii < M * N * W; ii++)
    {
        if (h.s[ii] < tol * maxs)
            donemask[ii] = 5;
    }

    LTFAT_NAME(borderstoheap)(&h, M, N, donemask);

    /* Create a struct holding inputs */
    hit = (struct LTFAT_NAME(heapinttask))
    {
        M, N, tgradw, fgradw, donemask
    };

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptorad)(tgrad, fgrad, a, M, L, tgradw, fgradw);

    while (1)
    {
        /* Inner loop processing all connected coefficients */
        while (h.heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(&h);
            /* Spread the current phase value to 4 direct neighbors */
            LTFAT_NAME(trapezheap)(&h, &hit, w, phase);
        }

        if (!LTFAT_NAME_REAL(findmaxinarraywrtmask)(h.s, donemask, M * N, &maxs, &Imax))
            break;

        /* Put maximal element onto the heap and mark that it is done. */
        LTFAT_NAME(heap_insert)(&h, Imax);
        donemask[Imax] = 6;
    }

    LTFAT_SAFEFREEALL(tgradw, fgradw, donemask, h.h, h.s);
}

void
LTFAT_NAME(borderstoheap)(struct LTFAT_NAME(heap)* h,
                          const ltfatInt M, const ltfatInt N,
                          int * donemask)
{
    for (ltfatInt w = 0; w < M * N ; w++)
    {
        // Is it a coefficient with known phase and is it big enough?
        // 5 is code of coefficients below tol
        if (donemask[w] && donemask[w] != 5)
        {
            // Is it a border coefficient?
            // i.e. is any of the 4 neighbors not reliable?
            if ( !donemask[NORTHFROMW(w, M, N)] ||
                    !donemask[ EASTFROMW(w, M, N)] ||
                    !donemask[SOUTHFROMW(w, M, N)] ||
                    !donemask[ WESTFROMW(w, M, N)] )
            {
                donemask[w] = 11; // Code of a good border coefficient
                LTFAT_NAME(heap_insert)(h, w);
            }
        }
    }
}



/*
 *  REAL-versions of the previous
 *
 *
 * */
void
LTFAT_NAME(borderstoheapreal)(struct LTFAT_NAME(heap)* h,
                              const ltfatInt M, const ltfatInt N,
                              int * donemask)
{
    ltfatInt M2 = M / 2 + 1;

    for (ltfatInt w = 0; w < M2 * N ; w++)
    {
        // Is it a coefficient with known phase and is it big enough?
        // 5 is code of coefficients below tol
        if (donemask[w] && donemask[w] != 5)
        {
            ltfatInt col = w / M2;
            ltfatInt row = w % M2;
            // Is it a border coefficient?
            // i.e. is any of the 4 neighbors not reliable?
            if ( ( row != M2 - 1   && !donemask[NORTHFROMW(w, M2, N)])  ||
                    ( col != N - 1 && !donemask[ EASTFROMW(w, M2, N)])  ||
                    ( row != 0     && !donemask[SOUTHFROMW(w, M2, N)])  ||
                    ( col != 0     && !donemask[ WESTFROMW(w, M2, N)]) )
            {
                donemask[w] = 11; // Code of a good border coefficient
                LTFAT_NAME(heap_insert)(h, w);
            }
        }
    }
}



void
LTFAT_NAME(gradsamptoradreal)(const LTFAT_REAL * tgrad, const LTFAT_REAL * fgrad,
                              ltfatInt a, ltfatInt M, ltfatInt L,
                              LTFAT_REAL * tgradw, LTFAT_REAL * fgradw)
{
    ltfatInt N = L / a;
    ltfatInt b = L / M;
    ltfatInt M2 = M / 2 + 1;
    LTFAT_REAL sampToRadConst = (LTFAT_REAL)( 2.0 * PI / L);

    for (ltfatInt jj = 0; jj < N; jj++)
    {
        for (ltfatInt ii = 0; ii < M2; ii++)
        {
            tgradw[ii + jj * M2] =    a * tgrad[ii + jj * M2] * sampToRadConst;
            fgradw[ii + jj * M2] =  - b * ( fgrad[ii + jj * M2] + jj * a ) * sampToRadConst;
            /* The following converts phase derivatives so that the result is time-invariant phase */
            /* tgradw[ii + jj * M] =    a * (tgrad[ii + jj * M] + ii*b) * sampToRadConst; */
            /* fgradw[ii + jj * M] =  - b * ( fgrad[ii + jj * M] ) * sampToRadConst; */
        }
    }
}



void LTFAT_NAME(trapezheapreal)(struct LTFAT_NAME(heap) *h,
                                const struct LTFAT_NAME(heapinttask) *heaptask,
                                const ltfatInt w,
                                LTFAT_REAL * phase)
{
    const ltfatInt M = heaptask->M;
    const ltfatInt N = heaptask->N;
    const LTFAT_REAL* tgradw = heaptask->tgrad;
    const LTFAT_REAL* fgradw = heaptask->fgrad;
    int* donemask = heaptask->donemask;
    ltfatInt w_E, w_W, w_N, w_S, row, col;

    /* North */
    w_N = NORTHFROMW(w, M, N);
    /* South */
    w_S = SOUTHFROMW(w, M, N);
    /* East */
    w_E = EASTFROMW(w, M, N);
    /* West */
    w_W = WESTFROMW(w, M, N);

    col = w / M;
    row = w % M;

    /* Try and put the four neighbours onto the heap.
     * Integration by trapezoidal rule */

    if (!donemask[w_N] && row != M - 1 )
    {
        phase[w_N] = phase[w] + (fgradw[w] + fgradw[w_N]) / 2;
        donemask[w_N] = 1;
        LTFAT_NAME(heap_insert)(h, w_N);
    }

    if (!donemask[w_S] && row != 0)
    {
        phase[w_S] = phase[w] - (fgradw[w] + fgradw[w_S]) / 2;
        donemask[w_S] = 2;
        LTFAT_NAME(heap_insert)(h, w_S);
    }

    if (!donemask[w_E] && col != N - 1)
    {
        phase[w_E] = phase[w] + (tgradw[w] + tgradw[w_E]) / 2;
        donemask[w_E] = 3;
        LTFAT_NAME(heap_insert)(h, w_E);
    }

    if (!donemask[w_W] && col != 0)
    {
        phase[w_W] = phase[w] - (tgradw[w] + tgradw[w_W]) / 2;
        donemask[w_W] = 4;
        LTFAT_NAME(heap_insert)(h, w_W);
    }
}


LTFAT_EXTERN
void LTFAT_NAME(heapintreal)(const LTFAT_REAL * s,
                             const LTFAT_REAL * tgrad,
                             const LTFAT_REAL * fgrad,
                             const ltfatInt a, const ltfatInt M,
                             const ltfatInt L, const ltfatInt W,
                             LTFAT_REAL tol, LTFAT_REAL * phase)
{

    /* Declarations */
    ltfatInt N, ii, Imax, M2;
    ltfatInt w;
    LTFAT_REAL maxs;
    int *donemask;
    struct LTFAT_NAME(heap) h;
    struct LTFAT_NAME(heapinttask) hit;

    M2 = M / 2 + 1;
    N = L / a;

    /* Main body */

    h.totalheapsize  = (ltfatInt)(M2 * log((LTFAT_REAL)M2));
    h.h              = ltfat_malloc(h.totalheapsize * sizeof * h.h);
    h.s              = (LTFAT_REAL*)s;
    h.heapsize       = 0;

    donemask = ltfat_malloc(M2 * N * W * sizeof * donemask);

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL *tgradw = ltfat_malloc(M2 * N * sizeof * tgradw);
    LTFAT_REAL *fgradw = ltfat_malloc(M2 * N * sizeof * fgradw);

    /* Set the phase to zero initially */
    memset(phase, 0, M2 * N * W * sizeof * phase);

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptoradreal)(tgrad, fgrad, a, M, L, tgradw, fgradw);

    /* We will start intergration from the biffest coefficient */
    LTFAT_NAME_REAL(findmaxinarray)(s, M2 * N, &maxs, &Imax);

    /* Mark all the small elements as done, they get zero phase.
     * (But should get random phase instead)
     * Code 5
     */
    for (ii = 0; ii < M2 * N * W; ii++)
    {
        if (s[ii] < tol * maxs)
            donemask[ii] = 5;
        else
            donemask[ii] = 0;
    }

    /* Create a struct holding inputs */
    hit = (struct LTFAT_NAME(heapinttask))
    {
        M2, N, tgradw, fgradw, donemask
    };

    /* Outer loop over islands of coefficients */
    do
    {
        /* Empty the heap */
        h.heapsize = 0;

        /* Put maximal element onto the heap and mark that it is done.
         * It gets a zero phase.  */
        LTFAT_NAME(heap_insert)(&h, Imax);
        donemask[Imax] = 6;

        /* Inner loop processing all connected coefficients */
        while (h.heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(&h);
            /* Spread the current phase value to 4 direct neighbors */
            LTFAT_NAME(trapezheapreal)(&h, &hit, w, phase);
        }

        /* Find the new maximal element */
        /* Break from the outer loop when there is none */
    }
    while ( LTFAT_NAME_REAL(findmaxinarraywrtmask)(s, donemask, M2 * N, &maxs, &Imax) );

    LTFAT_SAFEFREEALL(tgradw, fgradw, donemask, h.h);
}


LTFAT_EXTERN
void LTFAT_NAME(maskedheapintreal)(const LTFAT_COMPLEX * c,
                                   const LTFAT_REAL * tgrad,
                                   const LTFAT_REAL * fgrad,
                                   const int* mask,
                                   const ltfatInt a, const ltfatInt M,
                                   const ltfatInt L, const ltfatInt W,
                                   LTFAT_REAL tol, LTFAT_REAL * phase)
{

    /* Declarations */
    ltfatInt N, ii, Imax, M2;
    ltfatInt w;
    LTFAT_REAL maxs;
    int *donemask;
    struct LTFAT_NAME(heap) h;
    struct LTFAT_NAME(heapinttask) hit;

    M2 = M / 2 + 1;
    N = L / a;

    /* Main body */

    h.s              = ltfat_malloc(M2 * N * sizeof * h.s);
    h.totalheapsize  = (ltfatInt)(M2 * log((LTFAT_REAL)M2));
    h.h              = ltfat_malloc(h.totalheapsize * sizeof * h.h);
    h.heapsize       = 0;

    for (ltfatInt w = 0; w < M2 * N * W; w++)
    {
        h.s[w] = LTFAT_COMPLEXH(cabs)(c[w]);
    }

    donemask = ltfat_malloc(M2 * N * W * sizeof * donemask);

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL *tgradw = ltfat_malloc(M2 * N * sizeof * tgradw);
    LTFAT_REAL *fgradw = ltfat_malloc(M2 * N * sizeof * fgradw);

    /* Copy known phase */
    /* Code 12 */
    for (ltfatInt w = 0; w < M2 * N * W; w++)
    {
        if (mask[w])
        {
            phase[w]    = LTFAT_COMPLEXH(carg)(c[w]);
            donemask[w] = 12; /* Code of known phase */
        }
        else
        {
            phase[w]    = 0;
            donemask[w] = 0; /* Mask of unknown phase */
        }
    }


    /* We will start intergration from the biffest coefficient */
    LTFAT_NAME_REAL(findmaxinarray)(h.s, M2 * N, &maxs, &Imax);

    /* Mark all the small elements as done, they get zero phase.
     * (But should get random phase instead)
     * Code 5
     */
    for (ii = 0; ii < M2 * N * W; ii++)
    {
        if (h.s[ii] < tol * maxs)
            donemask[ii] = 5;
    }

    LTFAT_NAME(borderstoheapreal)(&h, M, N, donemask);

    /* Create a struct holding inputs */
    hit = (struct LTFAT_NAME(heapinttask))
    {
        M2, N, tgradw, fgradw, donemask
    };

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptoradreal)(tgrad, fgrad, a, M, L, tgradw, fgradw);

    while (1)
    {
        /* Inner loop processing all connected coefficients */
        while (h.heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(&h);
            /* Spread the current phase value to 4 direct neighbors */
            LTFAT_NAME(trapezheapreal)(&h, &hit, w, phase);
        }

        if (!LTFAT_NAME_REAL(findmaxinarraywrtmask)(h.s, donemask, M2 * N, &maxs, &Imax))
            break;

        /* Put maximal element onto the heap and mark that it is done. */
        LTFAT_NAME(heap_insert)(&h, Imax);
        donemask[Imax] = 6;
    }


    LTFAT_SAFEFREEALL(tgradw, fgradw, donemask, h.h, h.s);
}

#undef NORTHFROMW
#undef SOUTHFROMW
#undef WESTFROMW
#undef EASTFROMW


