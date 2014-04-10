#ifndef LTFAT_H
#define LTFAT_H 1
#include "config.h"

#ifdef LTFAT_COMPAT32
typedef int       ltfatInt;
#else
typedef ptrdiff_t ltfatInt;
#endif /* defined(LTFAT_COMPAT32) */

typedef enum
{
   ltfatUnspecErr = 1,
   ltfatNoErr = 0
} ltfatStatus;

typedef enum
{
   FREQINV = 0,
   TIMEINV = 1
} dgt_phasetype;

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

typedef enum
{
    CZT_NEXTFASTFFT,
    CZT_NEXTPOW2
} czt_ffthint;


/* BEGIN_C_DECLS */

#ifdef __cplusplus
extern "C"
{
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
#undef LTFAT_COMPLEXH

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


#ifdef __cplusplus
}  // extern "C"
#endif
/* END_C_DECLS */

#endif /* !LTFAT_H */
