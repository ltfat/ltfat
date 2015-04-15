#ifndef _LTFAT_TYPECONSTANT
#define _LTFAT_TYPECONSTANT
#include "config.h"



/* -------- Define routines that do not change between single/double-- */
LTFAT_EXTERN_TOO
ltfatInt gcd(const ltfatInt a, const ltfatInt b, ltfatInt *r, ltfatInt *s );

LTFAT_EXTERN_TOO
void* ltfat_malloc (size_t n);

LTFAT_EXTERN_TOO
void* ltfat_calloc (size_t nmemb, size_t size);

LTFAT_EXTERN_TOO
void* ltfat_realloc (void *ptr, size_t n);

LTFAT_EXTERN_TOO
void* ltfat_realloc_and_copy (void *ptr, size_t nold, size_t nnew);

LTFAT_EXTERN_TOO
void  ltfat_free(const void *ptr);

LTFAT_EXTERN_TOO
void  ltfat_safefree(const void *ptr);

LTFAT_EXTERN_TOO
void fftindex(const ltfatInt N, ltfatInt *indexout);

LTFAT_EXTERN_TOO
ltfatInt makelarger(const ltfatInt L, const ltfatInt K);

LTFAT_EXTERN_TOO
ltfatInt filterbank_td_size(const ltfatInt L,const ltfatInt a,
                            const ltfatInt gl,const ltfatInt offset,
                            const ltfatExtType ext);

LTFAT_EXTERN_TOO
ltfatInt imax(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN_TOO
ltfatInt imin(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN_TOO
ltfatInt lcm(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN_TOO
void gabimagepars(const ltfatInt Ls, const ltfatInt x, const ltfatInt y,
                  ltfatInt *a, ltfatInt *M, ltfatInt *L, ltfatInt *N, ltfatInt *Ngood);

LTFAT_EXTERN_TOO
ltfatInt wfacreal_size(const ltfatInt L, const ltfatInt a, const ltfatInt M);

LTFAT_EXTERN_TOO
ltfatInt nextPow2(const ltfatInt x);

LTFAT_EXTERN_TOO
ltfatInt nextfastfft(const ltfatInt x);

LTFAT_EXTERN_TOO
ltfatInt pow2(const ltfatInt x);

LTFAT_EXTERN_TOO
ltfatInt modPow2(const ltfatInt x,const ltfatInt pow2var);

LTFAT_EXTERN_TOO
ltfatInt ltfat_round(const double x);

LTFAT_EXTERN_TOO
ltfatInt positiverem(const ltfatInt a,const ltfatInt b);

LTFAT_EXTERN_TOO ltfatInt
rangelimit(const ltfatInt a, const ltfatInt amin, const ltfatInt amax);


// Custom headers are down here
#include "reassign_typeconstant.h"

#endif /* _LTFAT_TYPECONSTANT */
