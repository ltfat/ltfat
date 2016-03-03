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

LTFAT_EXTERN
void
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

    if (hit->do_real)
        LTFAT_NAME(borderstoheapreal)(hit->heap, hit->height, hit->N, hit->donemask);
    else
        LTFAT_NAME(borderstoheap)(hit->heap, hit->height, hit->N, hit->donemask);
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

#undef NORTHFROMW
#undef SOUTHFROMW
#undef WESTFROMW
#undef EASTFROMW
