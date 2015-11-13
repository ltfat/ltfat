/*This header file is suposed to be oncluded only from heapint.c */

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
    LTFAT_REAL *s;
};

struct LTFAT_NAME(heapinttask)
{
    ltfatInt M;
    ltfatInt N;
    const LTFAT_REAL *tgrad;
    const LTFAT_REAL *fgrad;
    int *donemask;
};


inline void
LTFAT_NAME(trapezheap)(struct LTFAT_NAME(heap) *heap,
                       const struct LTFAT_NAME(heapinttask) *heaptask,
                       const ltfatInt w, LTFAT_REAL* phase);

void
LTFAT_NAME(gradsamptorad)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                          ltfatInt a, ltfatInt M, ltfatInt L, dgt_phasetype phasetype,
                          LTFAT_REAL* tgradw, LTFAT_REAL* fgradw);

void
LTFAT_NAME(borderstoheap)(struct LTFAT_NAME(heap)* h,
                          const ltfatInt M, const ltfatInt N,
                          int * donemask);

inline void
LTFAT_NAME(trapezheapreal)(struct LTFAT_NAME(heap) *heap,
                           const struct LTFAT_NAME(heapinttask) *heaptask,
                           const ltfatInt w, LTFAT_REAL* phase);

void
LTFAT_NAME(gradsamptoradreal)(const LTFAT_REAL* tgrad, const LTFAT_REAL* fgrad,
                              ltfatInt a, ltfatInt M, ltfatInt L, dgt_phasetype phasetype,
                              LTFAT_REAL* tgradw, LTFAT_REAL* fgradw);

void
LTFAT_NAME(borderstoheapreal)(struct LTFAT_NAME(heap)* h,
                              const ltfatInt M, const ltfatInt N,
                              int * donemask);
