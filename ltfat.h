#ifndef LTFAT_H
#define LTFAT_H 1

#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "fftw3.h"
#include "cblas.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif /* defined(PI) */


#define LTFAT_MAKENAME(name,type,comp) name ## _ ## comp ## type
#define LTFAT_NAME_DOUBLE(name) LTFAT_MAKENAME(name,d,)
#define LTFAT_NAME_SINGLE(name) LTFAT_MAKENAME(name,s,)
#define LTFAT_NAME_COMPLEXDOUBLE(name) LTFAT_MAKENAME(name,d,c)
#define LTFAT_NAME_COMPLEXSINGLE(name) LTFAT_MAKENAME(name,s,c)

// "Vectorizes" a function call
#define LTFAT_APPLYFN(type,fn,...) do{ \
   const type list[] = {(const type)0,__VA_ARGS__}; \
   size_t len = sizeof(list)/sizeof(*list) - 1; \
   for(size_t ii=0;ii<len;ii++) \
      fn((const type)list[ii+1]); \
}while(0)

// Vectorized free
#define LTFAT_SAFEFREEALL(...) LTFAT_APPLYFN(void*,ltfat_safefree,__VA_ARGS__)

#ifdef LTFAT_COMPAT32
typedef int       ltfatInt;
#else
typedef ptrdiff_t ltfatInt;
#endif /* defined(PI) */

typedef enum
{
   ltfatUnspecErr = 1,
   ltfatNoErr = 0
} ltfatStatus;


typedef enum
{
    DCTI=FFTW_REDFT00, DCTIII=FFTW_REDFT01,
    DCTII=FFTW_REDFT10, DCTIV=FFTW_REDFT11
} dct_kind;


typedef enum
{
    DSTI=FFTW_RODFT00, DSTIII=FFTW_RODFT01,
    DSTII=FFTW_RODFT10, DSTIV=FFTW_RODFT11
} dst_kind;


/* BEGIN_C_DECLS */

#ifdef __cplusplus
extern "C"
{
#endif


/* Handle Windows DLL files */
/* defined by Makefile when compiling LTFAT */
#if defined(DLL_EXPORT_SYMBOLS) && (defined(_WIN32) || defined(__WIN32__))
#  define LTFAT_EXTERN extern __declspec(dllexport)
#else
#  define LTFAT_EXTERN
#endif


/* -------- Create the single precision routines headers ----- */

#ifndef LTFAT_DOUBLE
#ifndef LTFAT_SINGLE
#      define LTFAT_SINGLE_WASNOTDEFINED
#      define LTFAT_SINGLE
#endif

#   include "ltfat_types.h"
#   include "ltfat_typecomplexindependent.h"

#   define LTFAT_COMPLEXTYPE
#   include "ltfat_types.h"
#   include "ltfat_typecomplexindependent.h"
#   undef LTFAT_COMPLEXTYPE

#   include "ltfat_types.h"
#   include "ltfat_typeindependent.h"

#   ifdef LTFAT_SINGLE_WASNOTDEFINED
#      undef LTFAT_SINGLE
#      undef LTFAT_SINGLE_WASNOTDEFINED
#   endif
#endif


/* -------- Create the single precision routines headers ----- */
#ifndef LTFAT_SINGLE
#ifndef LTFAT_DOUBLE
#   define LTFAT_DOUBLE_WASNOTDEFINED
#   define LTFAT_DOUBLE
#endif

#include "ltfat_types.h"
#include "ltfat_typecomplexindependent.h"

#define LTFAT_COMPLEXTYPE
#include "ltfat_types.h"
#include "ltfat_typecomplexindependent.h"
#undef LTFAT_COMPLEXTYPE

#include "ltfat_types.h"
#include "ltfat_typeindependent.h"

#   ifdef LTFAT_DOUBLE_WASNOTDEFINED
#      undef LTFAT_DOUBLE
#      undef LTFAT_DOUBLE_WASNOTDEFINED
#   endif
#endif


// Undef all
#undef LTFAT_COMPLEX
#undef LTFAT_REAL
#undef LTFAT_TYPE
#undef LTFAT_NAME
#undef LTFAT_NAME_REAL
#undef LTFAT_NAME_COMPLEX
#undef LTFAT_FFTW
#undef LTFAT_MX_CLASSID
#undef LTFAT_MX_COMPLEXITY
#undef LTFAT_COMPLEXH_NAME


#ifdef LTFAT_DOUBLE
#define LTFAT_EXTERN_TOO LTFAT_EXTERN
#else
#define LTFAT_EXTERN_TOO
#endif


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
ltfatInt filterbank_td_size(ltfatInt L, ltfatInt a, ltfatInt gl, ltfatInt offset, ltfatExtType ext);

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
ltfatInt nextPow2(ltfatInt x);

LTFAT_EXTERN_TOO
ltfatInt nextfastfft(ltfatInt x);

LTFAT_EXTERN_TOO
ltfatInt pow2(ltfatInt x);

LTFAT_EXTERN_TOO
ltfatInt modPow2(ltfatInt x,ltfatInt pow2var);


#ifdef __cplusplus
}  // extern "C"
#endif
/* END_C_DECLS */

#endif /* !LTFAT_H */
