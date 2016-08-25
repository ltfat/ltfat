#ifndef _LTFAT_H
#define _LTFAT_H 1
//#include "ltfat/config.h"
#include "ltfat/complexcompat.h"

#ifndef NOSYSTEMHEADERS
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#endif


#ifdef LTFAT_COMPAT32
typedef int       ltfat_int;
#else
typedef ptrdiff_t ltfat_int;
#endif /* defined(LTFAT_COMPAT32) */


#include "ltfat/basicmacros.h"

/* BEGIN_C_DECLS */
#ifdef __cplusplus
extern "C"
{
#endif

#include "ltfat/errno.h"
#include "ltfat/version.h"

// First, include headers of type (single, double, or complex versions) inert functions
#include "ltfat/typeconstant.h"


/* -------- Create the single precision routines headers ----- */

#ifndef LTFAT_DOUBLE
#   ifndef LTFAT_SINGLE
#      define LTFAT_SINGLE_WASNOTDEFINED
#      define LTFAT_SINGLE
#   endif

#   include "ltfat/typecomplexindependent.h"

#   ifndef LTFAT_COMPLEXTYPE
#       define LTFAT_COMPLEXTYPE
#       include "ltfat/typecomplexindependent.h"
#       undef LTFAT_COMPLEXTYPE
#       include "ltfat/typeindependent.h"
#   else
#       undef LTFAT_COMPLEXTYPE
#       include "ltfat/typeindependent.h"
#       include "ltfat/typecomplexindependent.h"
#       define LTFAT_COMPLEXTYPE
#   endif

#   ifdef LTFAT_SINGLE_WASNOTDEFINED
#      undef LTFAT_SINGLE
#      undef LTFAT_SINGLE_WASNOTDEFINED
#   endif
#endif


/* -------- Create the single precision routines headers ----- */
#ifndef LTFAT_SINGLE
#   ifndef LTFAT_DOUBLE
#       define LTFAT_DOUBLE_WASNOTDEFINED
#       define LTFAT_DOUBLE
#   endif

#   include "ltfat/typecomplexindependent.h"

#   ifndef LTFAT_COMPLEXTYPE
#       define LTFAT_COMPLEXTYPE
#       include "ltfat/typecomplexindependent.h"
#       undef LTFAT_COMPLEXTYPE
#       include "ltfat/typeindependent.h"
#   else
#       undef LTFAT_COMPLEXTYPE
#       include "ltfat/typeindependent.h"
#       include "ltfat/typecomplexindependent.h"
#       define LTFAT_COMPLEXTYPE
#   endif

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
