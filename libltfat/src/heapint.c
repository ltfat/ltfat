#include "ltfat.h"
#include "ltfat_types.h"
#include "heapint.h"

#define NORTHFROMW(w,M,N) ((((w) + 1) % (M)) + (w) - (w) % (M))
#define SOUTHFROMW(w,M,N) (((w) - 1 + (M)) % (M) + (w) - (w) % (M))

#define EASTFROMW(w,M,N)  (((w) + (M)) % ((M) * (N)))
#define WESTFROMW(w,M,N)  (((w) - (M) + (M) * (N)) % ((M) * (N)))

struct LTFAT_NAME(heap)*
LTFAT_NAME(heap_init)(const ltfatInt initmaxsize, const LTFAT_REAL* s)
{
    struct LTFAT_NAME(heap)* h = ltfat_malloc(sizeof * h);

    h->totalheapsize  = initmaxsize;
    h->h              = ltfat_malloc(h->totalheapsize * sizeof * h->h);
    h->s              = s;
    h->heapsize       = 0;
    return h;
}

void
LTFAT_NAME(heap_done)(struct LTFAT_NAME(heap)* h)
{
    ltfat_free(h->h);
    ltfat_free(h);
}

void
LTFAT_NAME(heap_reset)(struct LTFAT_NAME(heap)* h, const LTFAT_REAL* news)
{
    h->s = news;
    h->heapsize = 0;
}

void
LTFAT_NAME(heap_grow)(struct LTFAT_NAME(heap)* h, int factor)
{
    h->totalheapsize *= factor;
    h->h = ltfat_realloc_and_copy(h->h,
                                  h->totalheapsize * sizeof * h->h / 2,
                                  h->totalheapsize * sizeof * h->h);
}



LTFAT_EXTERN void
LTFAT_NAME(heap_insert)(struct LTFAT_NAME(heap) *h, const ltfatInt key)
{
    ltfatInt pos, pos2, swap;

    /* Grow heap if necessary */
    if (h->totalheapsize == h->heapsize)
    {
        LTFAT_NAME(heap_grow)( h, 2);
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

LTFAT_EXTERN ltfatInt
LTFAT_NAME(heap_delete)(struct LTFAT_NAME(heap) *h)
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

LTFAT_EXTERN
struct LTFAT_NAME(heapinttask)*
LTFAT_NAME(heapinttask_init)(const ltfatInt height, const ltfatInt N,
                             const ltfatInt initheapsize,
                             const LTFAT_REAL* s, int do_real)
{
    struct LTFAT_NAME(heapinttask)* hit = ltfat_malloc(sizeof * hit);
    hit->height = height;
    hit->N = N;
    hit->donemask = ltfat_malloc(height * N * sizeof * hit->donemask);
    hit->heap = LTFAT_NAME(heap_init)(initheapsize, s);
    hit->do_real = do_real;

    if (do_real)
        hit->intfun = LTFAT_NAME(trapezheapreal);
    else
        hit->intfun = LTFAT_NAME(trapezheap);

    return hit;
}

LTFAT_EXTERN
void
LTFAT_NAME(heapinttask_done)(struct LTFAT_NAME(heapinttask)* hit)
{
    if (hit->heap)
        LTFAT_NAME(heap_done)(hit->heap);

    ltfat_free(hit->donemask);
    ltfat_free(hit);
}

LTFAT_EXTERN
void
LTFAT_NAME(heapinttask_resetmax)(struct LTFAT_NAME(heapinttask)* hit,
                                 const LTFAT_REAL* news,
                                 const LTFAT_REAL tol)
{
    ltfatInt Imax;
    LTFAT_REAL maxs;

    LTFAT_NAME(heap_reset)(hit->heap, news);

    // Find the biggest coefficient
    LTFAT_NAME_REAL(findmaxinarray)(news,  hit->height * hit->N , &maxs, &Imax);

    /* Mark all the small elements as done, they get zero phase.
     * Code 5
     */
    for (ltfatInt ii = 0; ii < hit->height * hit->N; ii++)
    {
        if (news[ii] <= tol * maxs)
            hit->donemask[ii] = 5;
        else
            hit->donemask[ii] = 0;
    }

    LTFAT_NAME(heap_insert)(hit->heap, Imax);
    hit->donemask[Imax] = 6;
}

LTFAT_EXTERN void
LTFAT_NAME(heapinttask_resetmask)(struct LTFAT_NAME(heapinttask)* hit,
                                  const int* mask,
                                  const LTFAT_REAL* news,
                                  const LTFAT_REAL tol,
                                  const int do_log)
{
    ltfatInt Imax;
    LTFAT_REAL maxs;

    LTFAT_NAME(heap_reset)(hit->heap, news);

    /* Copy known phase */
    /* Code 12 */
    for (ltfatInt w = 0; w < hit->height * hit->N; w++)
    {
        if (mask[w])
        {
            hit->donemask[w] = 12; /* Code of known phase */
        }
        else
        {
            /* phase[w]    = 0; */
            hit->donemask[w] = 0; /* Mask of unknown phase */
        }
    }

    /* We will start intergration from the biggest coefficient */
    LTFAT_NAME_REAL(findmaxinarray)(news, hit->height * hit->N, &maxs, &Imax);

    /* Mark all the small elements as done, they get zero phase.
     * (But should get random phase instead)
     * Code 5
     */

    if (do_log)
    {
        for (ltfatInt ii = 0; ii < hit->height * hit->N; ii++)
            if (news[ii] <= tol + maxs)
                hit->donemask[ii] = 5;
    }
    else
    {
        for (ltfatInt ii = 0; ii < hit->height * hit->N; ii++)
            if (news[ii] <= tol * maxs)
                hit->donemask[ii] = 5;
    }


}

void LTFAT_NAME(trapezheap)(const struct LTFAT_NAME(heapinttask) *hit,
                            const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                            const ltfatInt w,
                            LTFAT_REAL* phase)
{
    const ltfatInt M = hit->height;
    const ltfatInt N = hit->N;
    struct LTFAT_NAME(heap)* h = hit->heap;
    int* donemask = hit->donemask;
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


void
LTFAT_NAME(gradsamptorad)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                          ltfatInt a, ltfatInt M, ltfatInt L, ltfatInt W,
                          dgt_phasetype phasetype, int do_real,
                          LTFAT_REAL* tgradw, LTFAT_REAL* fgradw)
{
    ltfatInt N = L / a;
    LTFAT_REAL b = ((LTFAT_REAL) L) / M;
    LTFAT_REAL sampToRadConst = (LTFAT_REAL)( 2.0 * M_PI / L);

    ltfatInt height = do_real ? M / 2 + 1 : M;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* tgradchan = tgrad + w * height * N;
        const LTFAT_REAL* fgradchan = fgrad + w * height * N;
        LTFAT_REAL* tgradwchan = tgradw + w * height * N;
        LTFAT_REAL* fgradwchan = fgradw + w * height * N;

        for (ltfatInt n = 0; n < N; n++)
        {
            for (ltfatInt m = 0; m < height; m++)
            {
                if (phasetype == FREQINV)
                {
                    tgradwchan[m + n * height] =    a * tgradchan[m + n * height] * sampToRadConst;
                    fgradwchan[m + n * height] =  - b * ( fgradchan[m + n * height] + n * a ) *
                                                  sampToRadConst;
                }
                else if (phasetype == TIMEINV)
                {
                    tgradwchan[m + n * height] =    a * (tgradchan[m + n * height] + m * b) *
                                                    sampToRadConst;
                    fgradwchan[m + n * height] =  - b * ( fgradchan[m + n * height] ) *
                                                  sampToRadConst;
                }
            }
        }
    }
}

LTFAT_EXTERN
void LTFAT_NAME(heapint)(const LTFAT_REAL* s,
                         const LTFAT_REAL* tgradw,
                         const LTFAT_REAL* fgradw,
                         const ltfatInt a, const ltfatInt M,
                         const ltfatInt L, const ltfatInt W,
                         LTFAT_REAL tol,  LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask)* hit;

    // Width of s
    ltfatInt N = L / a;

    /* Set the phase to zero initially */
    memset(phase, 0, M * N * W * sizeof * phase);

    // Init plan
    hit = LTFAT_NAME(heapinttask_init)( M, N, M * log((double)M) , s, 0);

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * M * N;
        const LTFAT_REAL* tgradwchan = tgradw + w * M * N;
        const LTFAT_REAL* fgradwchan = fgradw + w * M * N;
        LTFAT_REAL* phasechan = phase + w * M * N;

        LTFAT_NAME(heapinttask_resetmax)(hit, schan, tol);

        LTFAT_NAME(heapint_execute)(hit, tgradwchan, fgradwchan, phasechan);
    }

    LTFAT_NAME(heapinttask_done)(hit);
}

LTFAT_EXTERN
void LTFAT_NAME(maskedheapint)(const LTFAT_REAL* s,
                               const LTFAT_REAL* tgradw,
                               const LTFAT_REAL* fgradw,
                               const int* mask,
                               const ltfatInt a, const ltfatInt M,
                               const ltfatInt L, const ltfatInt W,
                               LTFAT_REAL tol,
                               LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask)* hit;

    ltfatInt N = L / a;

    /* Main body */
    hit = LTFAT_NAME(heapinttask_init)( M, N, M * log((double)M) , s, 0);

    // Set all phases outside of the mask to zeros, do not modify the rest
    for (ltfatInt ii = 0; ii < M * N * W; ii++)
        if (!mask[ii])
            phase[ii] = 0;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * M * N;
        const LTFAT_REAL* tgradwchan = tgradw + w * M * N;
        const LTFAT_REAL* fgradwchan = fgradw + w * M * N;
        const int* maskchan = mask + w * M * N;
        LTFAT_REAL* phasechan = phase + w * M * N;

        // Empty heap and fill it with the border coefficients from the mask
        LTFAT_NAME(heapinttask_resetmask)(hit, maskchan, schan, tol, 0);
        LTFAT_NAME(borderstoheap)(hit->heap, hit->height, hit->N, hit->donemask);


        LTFAT_NAME(heapint_execute)(hit, tgradwchan, fgradwchan, phasechan);
    }

    LTFAT_NAME(heapinttask_done)(hit);

}

void
LTFAT_NAME(borderstoheap)(struct LTFAT_NAME(heap)* h,
                          const ltfatInt height, const ltfatInt N,
                          int* donemask)
{
    for (ltfatInt w = 0; w < height * N ; w++)
    {
        // Is it a coefficient with known phase and is it big enough?
        // 5 is code of coefficients below tol
        if (donemask[w] && donemask[w] != 5)
        {
            // Is it a border coefficient?
            // i.e. is any of the 4 neighbors not reliable?
            if ( !donemask[NORTHFROMW(w, height, N)] ||
                 !donemask[ EASTFROMW(w, height, N)] ||
                 !donemask[SOUTHFROMW(w, height, N)] ||
                 !donemask[ WESTFROMW(w, height, N)] )
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
                              const ltfatInt height, const ltfatInt N,
                              int* donemask)
{

    for (ltfatInt w = 0; w < height * N ; w++)
    {
        // Is it a coefficient with known phase and is it big enough?
        if (donemask[w] == 12)
        {
            ltfatInt col = w / height;
            ltfatInt row = w % height;

            // Is it a border coefficient?
            // i.e. is any of the 4 neighbors not reliable?
            if ( ( row != height - 1   && !donemask[NORTHFROMW(w, height, N)])  ||
                 ( col != N - 1    && !donemask[ EASTFROMW(w, height, N)])  ||
                 ( row != 0        && !donemask[SOUTHFROMW(w, height, N)])  ||
                 ( col != 0        && !donemask[ WESTFROMW(w, height, N)]) )
            {
                donemask[w] = 11; // Code of a good border coefficient
                LTFAT_NAME(heap_insert)(h, w);
            }
        }
    }
}

void
LTFAT_NAME(borderstoheapneighs)(struct LTFAT_NAME(heap)* h,
                                const ltfatInt Nsum, const ltfatInt neighs[], int* donemask)
{

    for (ltfatInt n = 0; n < Nsum; n++)
    {
        if (donemask[n] && donemask[n] != 5)
        {
            const ltfatInt* wneigh = neighs + 6 * n;

            for (ltfatInt ii = 0; ii < 6; ii++)
            {
                if ( !donemask[wneigh[ii]] )
                {
                    donemask[n] = 11; // Code of a good border coefficient
                    LTFAT_NAME(heap_insert)(h, n);
                    break;
                }
                // Do nothing if none of the neighbors is unknown
            }
        }
    }

}

void LTFAT_NAME(trapezheapreal)(const struct LTFAT_NAME(heapinttask) *hit,
                                const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                                const ltfatInt w,
                                LTFAT_REAL* phase)
{
    const ltfatInt M2 = hit->height;
    const ltfatInt N = hit->N;
    int* donemask = hit->donemask;
    struct LTFAT_NAME(heap) *h = hit->heap;
    ltfatInt w_E, w_W, w_N, w_S, row, col;

    /* North */
    w_N = NORTHFROMW(w, M2, N);
    /* South */
    w_S = SOUTHFROMW(w, M2, N);
    /* East */
    w_E = EASTFROMW(w, M2, N);
    /* West */
    w_W = WESTFROMW(w, M2, N);

    col = w / M2;
    row = w % M2;

    /* Try and put the four neighbours onto the heap.
     * Integration by trapezoidal rule */

    if (!donemask[w_N] && row != M2 - 1 )
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
void LTFAT_NAME(heapint_execute)(struct LTFAT_NAME(heapinttask)* hit,
                                 const LTFAT_REAL* tgradw,
                                 const LTFAT_REAL* fgradw,
                                 LTFAT_REAL* phase)
{
    /* Declarations */
    ltfatInt Imax;
    ltfatInt w;
    LTFAT_REAL maxs;
    int* donemask = hit->donemask;
    struct LTFAT_NAME(heap)* h = hit->heap;

    while (1)
    {
        /* Inner loop processing all connected coefficients */
        while (h->heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(h);
            /* Spread the current phase value to 4 direct neighbors */
            (*hit->intfun)(hit, tgradw, fgradw, w, phase);
        }

        if (!LTFAT_NAME_REAL(findmaxinarraywrtmask)(h->s, donemask,
                hit->height * hit->N, &maxs, &Imax))
            break;

        /* Put maximal element onto the heap and mark that it is done. */
        LTFAT_NAME(heap_insert)(h, Imax);
        donemask[Imax] = 6;
    }
}

/*
 * tgradw and fgradw must be in radians and scaled such that the step is 1
 */
LTFAT_EXTERN
void LTFAT_NAME(heapintreal)(const LTFAT_REAL* s,
                             const LTFAT_REAL* tgradw,
                             const LTFAT_REAL* fgradw,
                             const ltfatInt a, const ltfatInt M,
                             const ltfatInt L, const ltfatInt W,
                             LTFAT_REAL tol, LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask)* hit;

    // Height of s
    ltfatInt M2 = M / 2 + 1;
    // Width of s
    ltfatInt N = L / a;

    /* Set the phase to zero initially */
    memset(phase, 0, M2 * N * W * sizeof * phase);

    // Init plan
    hit = LTFAT_NAME(heapinttask_init)( M2, N, M2 * log((double)M2), s, 1);

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * M2 * N;
        const LTFAT_REAL* tgradwchan = tgradw + w * M2 * N;
        const LTFAT_REAL* fgradwchan = fgradw + w * M2 * N;
        LTFAT_REAL* phasechan = phase + w * M2 * N;

        // empty heap and add max element to it
        LTFAT_NAME(heapinttask_resetmax)(hit, schan, tol);

        LTFAT_NAME(heapint_execute)(hit, tgradwchan, fgradwchan, phasechan);
    }

    LTFAT_NAME(heapinttask_done)(hit);
}


LTFAT_EXTERN
void LTFAT_NAME(maskedheapintreal)(const LTFAT_REAL* s,
                                   const LTFAT_REAL* tgradw,
                                   const LTFAT_REAL* fgradw,
                                   const int* mask,
                                   const ltfatInt a, const ltfatInt M,
                                   const ltfatInt L, const ltfatInt W,
                                   LTFAT_REAL tol, LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask)* hit;

    ltfatInt M2 = M / 2 + 1;
    ltfatInt N = L / a;

    // Initialize plan
    hit = LTFAT_NAME(heapinttask_init)( M2, N, M2 * log((double) M2), s, 1);

    // Set all phases outside of the mask to zeros, do not modify the rest
    for (ltfatInt ii = 0; ii < M2 * N * W; ii++)
        if (!mask[ii])
            phase[ii] = 0;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * M2 * N;
        const LTFAT_REAL* tgradwchan = tgradw + w * M2 * N;
        const LTFAT_REAL* fgradwchan = fgradw + w * M2 * N;
        const int* maskchan = mask + w * M2 * N;
        LTFAT_REAL* phasechan = phase + w * M2 * N;

        // Empty heap and fill it with the border coefficients from the mask
        LTFAT_NAME(heapinttask_resetmask)(hit, maskchan, schan, tol, 0);
        LTFAT_NAME(borderstoheapreal)(hit->heap, hit->height, hit->N, hit->donemask);

        LTFAT_NAME(heapint_execute)(hit, tgradwchan, fgradwchan, phasechan);
    }

    LTFAT_NAME(heapinttask_done)(hit);
}

/*
 *  The _relgrad versions are just wrappers.
 *  They convert the relative phase gradients in samples to
 *  absolute phase gradinets in radians.
 * */

LTFAT_EXTERN void
LTFAT_NAME(maskedheapint_relgrad)(const LTFAT_REAL* s,
                                  const LTFAT_REAL* tgrad,
                                  const LTFAT_REAL* fgrad,
                                  const int* mask,
                                  const ltfatInt a, const ltfatInt M,
                                  const ltfatInt L, const ltfatInt W,
                                  const LTFAT_REAL tol, dgt_phasetype phasetype,
                                  LTFAT_REAL* phase)
{
    ltfatInt N = L / a;

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(M * N * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(M * N * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptorad)(tgrad, fgrad, a, M, L, W, phasetype, 0, tgradw,
                              fgradw);

    LTFAT_NAME(maskedheapint)(s, tgradw, fgradw, mask, a, M, L, W, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

LTFAT_EXTERN void
LTFAT_NAME(heapint_relgrad)(const LTFAT_REAL* s,
                            const LTFAT_REAL* tgrad,
                            const LTFAT_REAL* fgrad,
                            const ltfatInt a, const ltfatInt M,
                            const ltfatInt L, const ltfatInt W,
                            const LTFAT_REAL tol, dgt_phasetype phasetype,
                            LTFAT_REAL* phase)
{
    ltfatInt N = L / a;

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(M * N * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(M * N * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptorad)(tgrad, fgrad, a, M, L, W, phasetype, 0, tgradw,
                              fgradw);

    LTFAT_NAME(heapint)(s, tgradw, fgradw, a, M, L, W, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

LTFAT_EXTERN void
LTFAT_NAME(maskedheapintreal_relgrad)(const LTFAT_REAL* s,
                                      const LTFAT_REAL* tgrad,
                                      const LTFAT_REAL* fgrad,
                                      const int* mask,
                                      const ltfatInt a, const ltfatInt M,
                                      const ltfatInt L, const ltfatInt W,
                                      LTFAT_REAL tol, dgt_phasetype phasetype, LTFAT_REAL* phase)
{
    ltfatInt M2 = M / 2 + 1;
    ltfatInt N = L / a;

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(M2 * N * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(M2 * N * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptorad)(tgrad, fgrad, a, M, L, W, phasetype, 1, tgradw,
                              fgradw);

    LTFAT_NAME(maskedheapintreal)(s, tgradw, fgradw, mask, a, M, L, W, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);

}

LTFAT_EXTERN
void LTFAT_NAME(heapintreal_relgrad)(const LTFAT_REAL* s,
                                     const LTFAT_REAL* tgrad,
                                     const LTFAT_REAL* fgrad,
                                     const ltfatInt a, const ltfatInt M,
                                     const ltfatInt L, const ltfatInt W,
                                     LTFAT_REAL tol, dgt_phasetype phasetype, LTFAT_REAL* phase)
{
    ltfatInt M2 = M / 2 + 1;
    ltfatInt N = L / a;

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(M2 * N * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(M2 * N * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in both
     * directions */
    LTFAT_NAME(gradsamptorad)(tgrad, fgrad, a, M, L, W, phasetype, 1, tgradw,
                              fgradw);

    LTFAT_NAME(heapintreal)(s, tgradw, fgradw, a, M, L, W, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

/*--------------------------------FILTERBANK HEAP INTEGRATION---------------------------------*/
LTFAT_EXTERN
struct LTFAT_NAME(heapinttask_ufb)*
LTFAT_NAME(heapinttask_init_ufb)(const ltfatInt height, const ltfatInt N,
                                 const ltfatInt initheapsize,
                                 const LTFAT_REAL* s, int do_real)
{
    struct LTFAT_NAME(heapinttask_ufb)* fbhit = ltfat_malloc(sizeof * fbhit);
    fbhit->hit = LTFAT_NAME(heapinttask_init)( height, N, initheapsize, s, do_real);

    if (do_real)
        fbhit->intfun = LTFAT_NAME(trapezheapreal_ufb);
    else
        fbhit->intfun = LTFAT_NAME(trapezheap_ufb);

    return fbhit;
}

/* Execute the Heap Integration */
LTFAT_EXTERN
void LTFAT_NAME(heapint_execute_ufb)(struct LTFAT_NAME(heapinttask_ufb)* fbhit,
                                     const LTFAT_REAL* tgradw,
                                     const LTFAT_REAL* fgradw,
                                     const LTFAT_REAL* cfreq,
                                     LTFAT_REAL* phase)
{
    /* Declarations */
    ltfatInt Imax;
    ltfatInt w;
    LTFAT_REAL maxs;
    int* donemask = fbhit->hit->donemask;
    struct LTFAT_NAME(heap)* h = fbhit->hit->heap;

    while (1)
    {
        /* Inner loop processing all connected coefficients */
        while (h->heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(h);
            /* Spread the current phase value to 4 direct neighbors */
            (*fbhit->intfun)(fbhit->hit, tgradw, fgradw, cfreq, w, phase);
        }

        if (!LTFAT_NAME_REAL(findmaxinarraywrtmask)(h->s, donemask,
                fbhit->hit->height * fbhit->hit->N, &maxs, &Imax))
            break;

        /* Put maximal element onto the heap and mark that it is done. */
        LTFAT_NAME(heap_insert)(h, Imax);
        donemask[Imax] = 6;
    }
}

/* Trapezoidal integration rule, filterbank case with full frequency range */
void LTFAT_NAME(trapezheap_ufb)(const struct LTFAT_NAME(heapinttask) *hit,
                                const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                                const LTFAT_REAL* cfreq,
                                const ltfatInt w,
                                LTFAT_REAL* phase)
{
    const ltfatInt M = hit->height;
    const ltfatInt N = hit->N;
    struct LTFAT_NAME(heap)* h = hit->heap;
    int* donemask = hit->donemask;
    ltfatInt w_E, w_W, w_N, w_S, col, row;

    /* Try and put the four neighbours onto the heap.
     * Integration by trapezoidal rule */
    /* When integrating across frequencies, the difference of the associated
     * center frequencies has to be taken into account*/

    col = w / N;
    row = w % N;

    /* Inside a channel */

    /* South -> Backwards time */
    w_S = SOUTHFROMW(w, N, M);

    if (!donemask[w_S] && row != 0)
    {
        phase[w_S] = phase[w] - (tgradw[w] + tgradw[w_S]) / 2;
        donemask[w_S] = 3;
        LTFAT_NAME(heap_insert)(h, w_S);
    }

    /* North -> Forwards time */
    w_N = NORTHFROMW(w, N, M);

    if (!donemask[w_N] && row != N-1)
    {
        phase[w_N] = phase[w] + (tgradw[w] + tgradw[w_N]) / 2;
        donemask[w_N] = 4;
        LTFAT_NAME(heap_insert)(h, w_N);
    }

    /* Across channels */

    /* West -> Lower frequency */
    w_W = WESTFROMW(w, N, M);

    if (!donemask[w_W])
    {
        LTFAT_REAL step = cfreq[w_W / N] - cfreq[col];
        if (step > 0)
            step -= 2;

        phase[w_W] = phase[w] + step * (fgradw[w] + fgradw[w_W]) / 2;
        donemask[w_W] = 1;
        LTFAT_NAME(heap_insert)(h, w_W);
    }

    /* East -> Higher frequency */
    w_E = EASTFROMW(w, N, M);

    if (!donemask[w_E])
    {
        LTFAT_REAL step = cfreq[w_E / N] - cfreq[col];
        if (step < 0)
            step += 2;

        phase[w_E] = phase[w] + step * (fgradw[w] + fgradw[w_E]) / 2;
        donemask[w_E] = 2;
        LTFAT_NAME(heap_insert)(h, w_E);
    }    
}

/* Trapezoidal integration rule, filterbank case with partial frequency range -> Standard case*/
void LTFAT_NAME(trapezheapreal_ufb)(const struct LTFAT_NAME(heapinttask) *hit,
                                    const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                                    const LTFAT_REAL* cfreq,
                                    const ltfatInt w,
                                    LTFAT_REAL* phase)
{
    const ltfatInt M = hit->height;
    const ltfatInt N = hit->N;
    int* donemask = hit->donemask;
    struct LTFAT_NAME(heap) *h = hit->heap;
    ltfatInt w_E, w_W, w_N, w_S, row, col;    

    /* Try and put the four neighbours onto the heap.
     * Integration by trapezoidal rule */
    /* When integrating across frequencies, the difference of the associated
     * center frequencies has to be taken into account*/

    col = w / N;
    row = w % N;

    /* Inside a channel */

    /* South -> Backwards time */
    w_S = SOUTHFROMW(w, N, M);

    if (!donemask[w_S] && row != 0)
    {
        phase[w_S] = phase[w] - (tgradw[w] + tgradw[w_S]) / 2;
        donemask[w_S] = 3;
        LTFAT_NAME(heap_insert)(h, w_S);
    }

    /* North -> Forwards time*/
    w_N = NORTHFROMW(w, N, M);

    if (!donemask[w_N] && row != N-1)
    {
        phase[w_N] = phase[w] + (tgradw[w] + tgradw[w_N]) / 2;
        donemask[w_N] = 4;
        LTFAT_NAME(heap_insert)(h, w_N);
    }

    /* Across channels */

    /* West -> Lower frequency */
    w_W = WESTFROMW(w, N, M);

    if (!donemask[w_W] && col != 0)
    {
        LTFAT_REAL step = cfreq[w_W / N] - cfreq[col];
        if (step > 0)
            step -= 2;

        phase[w_W] = phase[w] + step * (fgradw[w] + fgradw[w_W]) / 2;
        donemask[w_W] = 1;
        LTFAT_NAME(heap_insert)(h, w_W);
    }

    /* East -> Higher frequency */
    w_E = EASTFROMW(w, N, M);

    if (!donemask[w_E] && col != M-1)
    {
        LTFAT_REAL step = cfreq[w_E / N] - cfreq[col];
        if (step < 0)
            step += 2;

        phase[w_E] = phase[w] + step * (fgradw[w] + fgradw[w_E]) / 2;
        donemask[w_E] = 2;
        LTFAT_NAME(heap_insert)(h, w_E);
    }    

}

/* Conversion of tgrad and fgrad to correct convention and scaling */
void
LTFAT_NAME(gradsamptorad_ufb)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                              const LTFAT_REAL* cfreq,
                              ltfatInt a, ltfatInt M, ltfatInt L, ltfatInt W,
                              LTFAT_REAL* tgradw, LTFAT_REAL* fgradw)
{
    ltfatInt N = L / a;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* tgradchan = tgrad + w * M * N;
        const LTFAT_REAL* fgradchan = fgrad + w * M * N;
        LTFAT_REAL* tgradwchan = tgradw + w * M * N;
        LTFAT_REAL* fgradwchan = fgradw + w * M * N;

        for (ltfatInt m = 0; m < M; m++)
        {
            for (ltfatInt n = 0; n < N; n++)
            {
                /*In contrast to Gabor, tgrad is not in samples, but in ]-1,1]*/
                tgradwchan[n + m * N] =    a * (tgradchan[n + m * N] + cfreq[m]) *
                                           M_PI;
                /*In contrast to Gabor, fgrad has to be weighted by the channel difference
                *DURING the integration. However, cfreq ranges in ]-1,1], so fgrad is 			*only scaled by PI.*/
                fgradwchan[n + m * N] =  - ( fgradchan[n + m * N] ) * M_PI;
            }
        }
    }
}

/* Interfacing functions for the various cases*/

LTFAT_EXTERN
void LTFAT_NAME(ufilterbankheapint)(const LTFAT_REAL* s,
                             const LTFAT_REAL* tgradw,
                             const LTFAT_REAL* fgradw,
                             const LTFAT_REAL* cfreq,
                             const ltfatInt a, const ltfatInt M,
                             const ltfatInt L, const ltfatInt W,
			     const int do_real, const LTFAT_REAL tol,  
		             LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask_ufb)* fbhit;

    // Width of s
    ltfatInt N = L / a;

    /* Set the phase to zero initially */
    memset(phase, 0, M * N * W * sizeof * phase);

    // Init plan
    fbhit = LTFAT_NAME(heapinttask_init_ufb)( M, N, M * log((double)M) , s, do_real);

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * M * N;
        const LTFAT_REAL* tgradwchan = tgradw + w * M * N;
        const LTFAT_REAL* fgradwchan = fgradw + w * M * N;
        LTFAT_REAL* phasechan = phase + w * M * N;

        LTFAT_NAME(heapinttask_resetmax)(fbhit->hit, schan, tol);

        LTFAT_NAME(heapint_execute_ufb)(fbhit, tgradwchan, fgradwchan, cfreq,
                                        phasechan);
    }

    LTFAT_NAME(heapinttask_done)(fbhit->hit);
    ltfat_free(fbhit);
}

LTFAT_EXTERN
void LTFAT_NAME(ufilterbankmaskedheapint)(const LTFAT_REAL* s,
                                   const LTFAT_REAL* tgradw,
                                   const LTFAT_REAL* fgradw,
                                   const LTFAT_REAL* cfreq,
                                   const int* mask,
                                   const ltfatInt a, const ltfatInt M,
                                   const ltfatInt L, const ltfatInt W,
				   const int do_real,
                                   const LTFAT_REAL tol,
                                   LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask_ufb)* fbhit;

    ltfatInt N = L / a;

    /* Main body */
    fbhit = LTFAT_NAME(heapinttask_init_ufb)( M, N, M * log((double)M) , s, do_real);

    // Set all phases outside of the mask to zeros, do not modify the rest
    for (ltfatInt ii = 0; ii < M * N * W; ii++)
        if (!mask[ii])
            phase[ii] = 0;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * M * N;
        const LTFAT_REAL* tgradwchan = tgradw + w * M * N;
        const LTFAT_REAL* fgradwchan = fgradw + w * M * N;
        const int* maskchan = mask + w * M * N;
        LTFAT_REAL* phasechan = phase + w * M * N;

        // Empty heap and fill it with the border coefficients from the mask
        LTFAT_NAME(heapinttask_resetmask)(fbhit->hit, maskchan, schan, tol, 0);
        LTFAT_NAME(borderstoheap)(fbhit->hit->heap, fbhit->hit->N, fbhit->hit->height,
                                  fbhit->hit->donemask);

        LTFAT_NAME(heapint_execute_ufb)(fbhit, tgradwchan, fgradwchan, cfreq,
                                        phasechan);
    }

    LTFAT_NAME(heapinttask_done)(fbhit->hit);
    ltfat_free(fbhit);

}

	/*  The _relgrad versions are just wrappers.
	 *  They convert the relative phase gradients in samples to
	 *  absolute phase gradinets in radians. */

LTFAT_EXTERN void
LTFAT_NAME(ufilterbankheapint_relgrad)(const LTFAT_REAL* s,
                                const LTFAT_REAL* tgrad,
                                const LTFAT_REAL* fgrad,
                                const LTFAT_REAL* cfreq,
                                const ltfatInt a, const ltfatInt M,
                                const ltfatInt L, const ltfatInt W,
			        const int do_real,
                                const LTFAT_REAL tol,
                                LTFAT_REAL* phase)
{
    ltfatInt N = L / a;

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(M * N * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(M * N * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in time
     * direction. The step in frequency direction is multiplied by the difference of the center
     * frequencies during integration.*/
    LTFAT_NAME(gradsamptorad_ufb)(tgrad, fgrad, cfreq, a, M, L, W, tgradw, fgradw);

    LTFAT_NAME(ufilterbankheapint)(s, tgradw, fgradw, cfreq, a, M, L, W, do_real, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

LTFAT_EXTERN void
LTFAT_NAME(ufilterbankmaskedheapint_relgrad)(const LTFAT_REAL* s,
                                      const LTFAT_REAL* tgrad,
                                      const LTFAT_REAL* fgrad,
                                      const LTFAT_REAL* cfreq,
                                      const int* mask,
                                      const ltfatInt a, const ltfatInt M,
                                      const ltfatInt L, const ltfatInt W,
				      const int do_real,
                                      const LTFAT_REAL tol,
                                      LTFAT_REAL* phase)
{
    ltfatInt N = L / a;

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(M * N * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(M * N * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in time
     * direction. The step in frequency direction is multiplied by the difference of the center
     * frequencies during integration.*/
    LTFAT_NAME(gradsamptorad_ufb)(tgrad, fgrad, cfreq, a, M, L, W, tgradw, fgradw);

    LTFAT_NAME(ufilterbankmaskedheapint)(s, tgradw, fgradw, cfreq, mask, a, M, L, W, do_real, tol,
                                  phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

/*--------------------------------GENERAL FILTER BANKS---------------------------------*/
LTFAT_EXTERN
struct LTFAT_NAME(heapinttask_fb)*
LTFAT_NAME(heapinttask_init_fb)(const ltfatInt height,
                                const ltfatInt initheapsize,
                                const LTFAT_REAL* s,
                                const ltfatInt* N,
                                const double* a,
                                const LTFAT_REAL* cfreq,
                                const ltfatInt* neigh,
                                const LTFAT_REAL* posInfo,
                                int do_real)
{
    struct LTFAT_NAME(heapinttask_fb)* fbhit = ltfat_malloc(sizeof * fbhit);
    fbhit->hit = LTFAT_NAME(heapinttask_init)( height, 1, initheapsize, s, do_real);

    fbhit->intfun = LTFAT_NAME(trapezheap_fb);

    fbhit->N = N;
    fbhit->a = a;
    fbhit->cfreq = cfreq;
    fbhit->neigh = neigh;
    fbhit->posInfo = posInfo;

    return fbhit;
}

/* Execute the Heap Integration */
LTFAT_EXTERN
void LTFAT_NAME(heapint_execute_fb)(struct LTFAT_NAME(heapinttask_fb)* fbhit,
                                    const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                                    LTFAT_REAL* phase)
{
    /* Declarations */
    ltfatInt Imax;
    ltfatInt w;
    LTFAT_REAL maxs;
    int* donemask = fbhit->hit->donemask;
    struct LTFAT_NAME(heap)* h = fbhit->hit->heap;
    
    while (1)
    {
        /* Inner loop processing all connected coefficients */
        while (h->heapsize)
        {
            /* Extract largest (first) element from heap and delete it. */
            w = LTFAT_NAME(heap_delete)(h);

            /* Spread the current phase value to 4 direct neighbors */
            (*fbhit->intfun)(fbhit, tgradw, fgradw, w, phase);
        }

        if (!LTFAT_NAME_REAL(findmaxinarraywrtmask)(h->s, donemask,
                fbhit->hit->height * fbhit->hit->N, &maxs, &Imax))
            break;
        /* Put maximal element onto the heap and mark that it is done. */
        LTFAT_NAME(heap_insert)(h, Imax);
        donemask[Imax] = 6;
    }
}

/* Trapeziodal rule for general filter banks */
void LTFAT_NAME(trapezheap_fb)(const struct LTFAT_NAME(heapinttask_fb) *fbhit,
                               const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,

                               const ltfatInt w, LTFAT_REAL* phase)
{
    struct LTFAT_NAME(heapinttask) * hit = fbhit->hit;
    struct LTFAT_NAME(heap)* h = hit->heap;
    int* donemask = hit->donemask;
    ltfatInt w_TMP;

    const ltfatInt* wneigh = fbhit->neigh + 6 * w;
    const LTFAT_REAL* posInfo = fbhit->posInfo;
    double* a = fbhit->a;
    const LTFAT_REAL* cfreq = fbhit->cfreq;

    /* Try and put all neighbors onto the heap, starting with neighbors in
     * the same channel, then next lower channel, finally next higher channel.
     * Integration by trapezoidal rule */
    /* When integrating across frequencies, the difference of the associated
     * center frequencies has to be taken into account*/

    /* Inside the channel */

    ltfatInt wchan = posInfo[2 * w];
    ltfatInt ii = 0;

    for (ii = 0; ii < 2; ++ii )
    {
        w_TMP = wneigh[ii];
        if (w_TMP >= 0 && !donemask[w_TMP])
        {
            phase[w_TMP] = phase[w] + a[wchan] * (w_TMP - w) *
                           (tgradw[w] + tgradw[w_TMP]) / 2;
            donemask[w_TMP] = 3;
            LTFAT_NAME(heap_insert)(h, w_TMP);
        }
    }

    /* Channel below */
    for (ii = 2; ii < 6; ++ii)
    {
        w_TMP = wneigh[ii];
        if (w_TMP >= 0 && !donemask[w_TMP])
        {
            phase[w_TMP] = phase[w] +
                           (posInfo[2 * w_TMP + 1] - posInfo[2 * w + 1] ) * (tgradw[w] + tgradw[w_TMP]) / 2
                           + (cfreq[(ltfatInt)posInfo[2 * w_TMP]] - cfreq[(ltfatInt)posInfo[2 * w]]) *
                           (fgradw[w] + fgradw[w_TMP]) / 2;
            donemask[w_TMP] = 3;
            LTFAT_NAME(heap_insert)(h, w_TMP);
        }
    }
}

/* Conversion of tgrad and fgrad to correct convention and scaling */
void
LTFAT_NAME(gradsamptorad_fb)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                             const LTFAT_REAL* cfreq,
                             const ltfatInt M,
                             const ltfatInt N[], const ltfatInt Nsum,
                             const ltfatInt W,
                             LTFAT_REAL* tgradw, LTFAT_REAL* fgradw)
{

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* tgradchan = tgrad + w * Nsum;
        const LTFAT_REAL* fgradchan = fgrad + w * Nsum;
        LTFAT_REAL* tgradwchan = tgradw + w * Nsum;
        LTFAT_REAL* fgradwchan = fgradw + w * Nsum;

        ltfatInt chanStart = 0;
        for (ltfatInt m = 0; m < M; m++)
        {
            for (ltfatInt n = 0; n < N[m]; n++)
            {                
                /*In contrast to Gabor, tgrad is not in samples, but in ]-1,1]*/
                tgradwchan[n + chanStart] =    (tgradchan[n + chanStart] +
                                                cfreq[m]) * M_PI;
                /*In contrast to Gabor, fgrad has to be weighted by the channel difference
                *DURING the integration. However, cfreq ranges in ]-1,1], so fgrad is
                *only scaled by PI.*/
                fgradwchan[n + chanStart] =  - ( fgradchan[n + chanStart] ) * M_PI;
            }
            chanStart += N[m];
        }
    }
}

/* Interfacing functions for the various cases*/

LTFAT_EXTERN
void LTFAT_NAME(filterbankheapint)(const LTFAT_REAL* s,
                            const LTFAT_REAL* tgradw,
                            const LTFAT_REAL* fgradw,
                            const ltfatInt neigh[],
                            const LTFAT_REAL posInfo[],
                            const LTFAT_REAL cfreq[],
                            const double a[], const ltfatInt M, const ltfatInt N[],
                            const ltfatInt Nsum, const ltfatInt W,
                            LTFAT_REAL tol,  LTFAT_REAL* phase)
{
    /* Declarations */
    struct LTFAT_NAME(heapinttask_fb)* fbhit;

    /* Set the phase to zero initially */
    memset(phase, 0, Nsum * W * sizeof * phase);

    // Init plan
    fbhit = LTFAT_NAME(heapinttask_init_fb)( Nsum, M * log((double)M) , s, N, a,
            cfreq, neigh, posInfo, 0);

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * Nsum;
        const LTFAT_REAL* tgradwchan = tgradw + w * Nsum;
        const LTFAT_REAL* fgradwchan = fgradw + w * Nsum;
        LTFAT_REAL* phasechan = phase + w * Nsum;

        LTFAT_NAME(heapinttask_resetmax)(fbhit->hit, schan, tol);

        LTFAT_NAME(heapint_execute_fb)(fbhit, tgradwchan, fgradwchan, phasechan);
    }

    LTFAT_NAME(heapinttask_done)(fbhit->hit);
    ltfat_free(fbhit);
}

LTFAT_EXTERN void
LTFAT_NAME(filterbankmaskedheapint)(const LTFAT_REAL* s,
                             const LTFAT_REAL* tgradw,
                             const LTFAT_REAL* fgradw,
                             const int* mask,
                             const ltfatInt neigh[],
                             const LTFAT_REAL posInfo[],
                             const LTFAT_REAL* cfreq,
                             const double* a,
                             const ltfatInt M,
                             const ltfatInt N[], const ltfatInt Nsum,
                             const ltfatInt W,
                             LTFAT_REAL tol,  LTFAT_REAL* phase)
{

    /* Declarations */
    struct LTFAT_NAME(heapinttask_fb)* fbhit;

    // Init plan
    fbhit = LTFAT_NAME(heapinttask_init_fb)( Nsum, M * log((double)M), s, N, a,
            cfreq, neigh, posInfo, 0);

    // Set all phases outside of the mask to zeros, do not modify the rest
    for (ltfatInt ii = 0; ii < W * Nsum; ii++)
        if (!mask[ii])
            phase[ii] = 0;

    for (ltfatInt w = 0; w < W; ++w)
    {
        const LTFAT_REAL* schan = s + w * Nsum;
        const LTFAT_REAL* tgradwchan = tgradw + w * Nsum;
        const LTFAT_REAL* fgradwchan = fgradw + w * Nsum;
        LTFAT_REAL* phasechan = phase + w * Nsum;
        const int* maskchan = mask + w * Nsum;

        LTFAT_NAME(heapinttask_resetmask)(fbhit->hit, maskchan, schan, tol, 0);
        LTFAT_NAME(borderstoheapneighs)(fbhit->hit->heap, Nsum, neigh,
                                        fbhit->hit->donemask);

        LTFAT_NAME(heapint_execute_fb)(fbhit, tgradwchan, fgradwchan, phasechan);
    }

    LTFAT_NAME(heapinttask_done)(fbhit->hit);
    ltfat_free(fbhit);
}

/*
 *  The _relgrad versions are just wrappers.
 *  They convert the relative phase gradients in samples to
 *  absolute phase gradinets in radians.
 * */

LTFAT_EXTERN void
LTFAT_NAME(filterbankheapint_relgrad)(const LTFAT_REAL* s,
                               const LTFAT_REAL* tgrad,
                               const LTFAT_REAL* fgrad,
                               const ltfatInt* neigh,
                               const LTFAT_REAL* posInfo,
                               const LTFAT_REAL* cfreq,
                               const double* a, const ltfatInt M,
                               const ltfatInt N[], const ltfatInt Nsum,
                               const ltfatInt W,
                               LTFAT_REAL tol,  LTFAT_REAL* phase)
{
    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(Nsum * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(Nsum * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in time
     * direction. The step in frequency direction is multiplied by the difference of the center
     * frequencies during integration.*/
    LTFAT_NAME(gradsamptorad_fb)(tgrad, fgrad, cfreq, M, N, Nsum, W,
                                 tgradw, fgradw);

    LTFAT_NAME(filterbankheapint)(s, tgradw, fgradw, neigh, posInfo, cfreq, a, M, N, Nsum,
                           W, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

LTFAT_EXTERN void
LTFAT_NAME(filterbankmaskedheapint_relgrad)(const LTFAT_REAL* s,
                                     const LTFAT_REAL* tgrad,
                                     const LTFAT_REAL* fgrad,
                                     const int* mask,
                                     const ltfatInt neigh[],
                                     const LTFAT_REAL posInfo[],
                                     const LTFAT_REAL* cfreq,
                                     const double* a,
                                     const ltfatInt M,
                                     const ltfatInt N[], const ltfatInt Nsum,
                                     const ltfatInt W,
                                     LTFAT_REAL tol,  LTFAT_REAL* phase)
{

    /* Allocate new arrays, we need to rescale the derivatives */
    LTFAT_REAL* tgradw = ltfat_malloc(Nsum * W * sizeof * tgradw);
    LTFAT_REAL* fgradw = ltfat_malloc(Nsum * W * sizeof * fgradw);

    /* Rescale the derivatives such that they are in radians and the step is 1 in time
     * direction. The step in frequency direction is multiplied by the difference of the center
     * frequencies during integration.*/
    LTFAT_NAME(gradsamptorad_fb)(tgrad, fgrad, cfreq, M, N, Nsum, W,
                                 tgradw, fgradw);

    LTFAT_NAME(filterbankmaskedheapint)(s, tgradw, fgradw, mask, neigh, posInfo, cfreq, a, M,
                                 N, Nsum, W, tol, phase);

    LTFAT_SAFEFREEALL(tgradw, fgradw);
}

#undef NORTHFROMW
#undef SOUTHFROMW
#undef WESTFROMW
#undef EASTFROMW
