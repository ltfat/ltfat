/*This header file is suposed to be oncluded only from heapint.c */

/*
 * Heap is a dynamic array h of heapsize elements which are ordered according
 * to s such that s[h[0]] is the maximum.
 *
 * */
#include "ltfat.h"
#include "ltfat/types.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct LTFAT_NAME(heap) LTFAT_NAME(heap);

LTFAT_NAME(heap)*
LTFAT_NAME(heap_init)(ltfat_int initmaxsize, const LTFAT_REAL* s);

void
LTFAT_NAME(heap_done)(LTFAT_NAME(heap)* h);

void
LTFAT_NAME(heap_grow)(LTFAT_NAME(heap)* h, int factor);

void
LTFAT_NAME(heap_reset)(LTFAT_NAME(heap)* h, const LTFAT_REAL* news);

LTFAT_API ltfat_int
LTFAT_NAME(heap_delete)(LTFAT_NAME(heap) *h);

LTFAT_API void
LTFAT_NAME(heap_insert)(LTFAT_NAME(heap) *h, ltfat_int key);

/*  */

inline void
LTFAT_NAME(trapezheap)(const LTFAT_NAME(heapinttask) *heaptask,
                       const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                       ltfat_int w, LTFAT_REAL* phase);

inline void
LTFAT_NAME(trapezheapreal)(const LTFAT_NAME(heapinttask) *heaptask,
                           const LTFAT_REAL* tgradw, const LTFAT_REAL* fgradw,
                           ltfat_int w, LTFAT_REAL* phase);

void
LTFAT_NAME(gradsamptorad)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                          ltfat_int a, ltfat_int M, ltfat_int L, ltfat_int W,
                          ltfat_phaseconvention phasetype, int do_real,
                          LTFAT_REAL* tgradw, LTFAT_REAL* fgradw);

void
LTFAT_NAME(borderstoheap)(LTFAT_NAME(heap)* h,
                          ltfat_int height, ltfat_int N,
                          int * donemask);

void
LTFAT_NAME(borderstoheapreal)(LTFAT_NAME(heap)* h,
                              ltfat_int height, ltfat_int N,
                              int * donemask);

#ifdef __cplusplus
}
#endif
