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

typedef enum
{
    PER,
    PERDEC,
    PPD,
    SYM,
    EVEN,
    SYMW,
    ASYM,
    ODD,
    ASYMW,
    SP0,
    ZPD,
    ZERO,
    VALID,
    BAD_TYPE
} ltfatExtType;


/* BEGIN_C_DECLS */

#ifdef __cplusplus
extern "C"
{
#endif

// First, include headers of type (single, double, or complex versions) inert functions
#include "ltfat_typeconstant.h"

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


#ifdef __cplusplus
}  // extern "C"
#endif
/* END_C_DECLS */

#endif /* !LTFAT_H */
