#ifndef _LTFAT_TYPECONSTANT
#define _LTFAT_TYPECONSTANT
#include "config.h"



/* -------- Define routines that do not change between single/double-- */
LTFAT_EXTERN
ltfatInt gcd(const ltfatInt a, const ltfatInt b, ltfatInt *r, ltfatInt *s );

LTFAT_EXTERN
void* ltfat_malloc (size_t n);

LTFAT_EXTERN
void* ltfat_calloc (size_t nmemb, size_t size);

LTFAT_EXTERN
void* ltfat_realloc (void *ptr, size_t n);

LTFAT_EXTERN
void* ltfat_realloc_and_copy (void *ptr, size_t nold, size_t nnew);

LTFAT_EXTERN
void  ltfat_free(const void *ptr);

LTFAT_EXTERN
void  ltfat_safefree(const void *ptr);

LTFAT_EXTERN
void fftindex(const ltfatInt N, ltfatInt *indexout);

LTFAT_EXTERN
ltfatInt makelarger(const ltfatInt L, const ltfatInt K);

LTFAT_EXTERN
ltfatInt filterbank_td_size(const ltfatInt L,const ltfatInt a,
                            const ltfatInt gl,const ltfatInt offset,
                            const ltfatExtType ext);

LTFAT_EXTERN
ltfatInt imax(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN
ltfatInt imin(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN
ltfatInt lcm(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN
void gabimagepars(const ltfatInt Ls, const ltfatInt x, const ltfatInt y,
                  ltfatInt *a, ltfatInt *M, ltfatInt *L, ltfatInt *N, ltfatInt *Ngood);

LTFAT_EXTERN
ltfatInt wfacreal_size(const ltfatInt L, const ltfatInt a, const ltfatInt M);

LTFAT_EXTERN
ltfatInt nextPow2(const ltfatInt x);

LTFAT_EXTERN
ltfatInt nextfastfft(const ltfatInt x);

LTFAT_EXTERN
ltfatInt pow2(const ltfatInt x);

LTFAT_EXTERN
ltfatInt modPow2(const ltfatInt x,const ltfatInt pow2var);

LTFAT_EXTERN
ltfatInt ltfat_round(const double x);

LTFAT_EXTERN
ltfatInt positiverem(const ltfatInt a,const ltfatInt b);

LTFAT_EXTERN ltfatInt
rangelimit(const ltfatInt a, const ltfatInt amin, const ltfatInt amax);


// Custom headers are down here
#include "reassign_typeconstant.h"

#endif /* _LTFAT_TYPECONSTANT */
