// maxtree

typedef struct LTFAT_NAME(maxtree) LTFAT_NAME(maxtree);

LTFAT_API int
LTFAT_NAME(maxtree_init)( ltfat_int L, ltfat_int Lstep, ltfat_int depth, LTFAT_NAME(maxtree)** p);

LTFAT_API int
LTFAT_NAME(maxtree_initwitharray)(
    ltfat_int L, ltfat_int depth,
    const LTFAT_REAL* inarray,
    LTFAT_NAME(maxtree)** p);

LTFAT_API int
LTFAT_NAME(maxtree_reset)(LTFAT_NAME(maxtree)* p, const LTFAT_REAL* inarray);

LTFAT_API int
LTFAT_NAME(maxtree_setdirty)(LTFAT_NAME(maxtree)* p, ltfat_int start, ltfat_int end);

LTFAT_API int
LTFAT_NAME(maxtree_getdirty)(LTFAT_NAME(maxtree)* p, ltfat_int* start, ltfat_int* end);

LTFAT_API int
LTFAT_NAME(maxtree_updatedirty)(LTFAT_NAME(maxtree)* p);

LTFAT_API int
LTFAT_NAME(maxtree_updaterange)(LTFAT_NAME(maxtree)* p, ltfat_int start, ltfat_int stop);

LTFAT_API int
LTFAT_NAME(maxtree_findmax)(LTFAT_NAME(maxtree)* p, LTFAT_REAL* max, ltfat_int* maxPos);

LTFAT_API int
LTFAT_NAME(maxtree_done)(LTFAT_NAME(maxtree)** p);
