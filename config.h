#ifndef CONFIG_H
#define CONFIG_H 1
/**
* Include files
* Constants
* Macros, not changing
*/
#ifdef __cplusplus
   // C++ complex header
   // fftw3.h will define:
   // typedef double fftw_complex[2]
   #include <complex>
#else
   // C99 complex header
   // fftw3.h will define:
   // typedef double _Complex fftw_complex
   #include <complex.h>
#endif
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
//#include <stdbool.h> // We do not use bool anyway
#include <math.h>
#include "fftw3.h"
#include "cblas.h"

#define HAVE_BLAS 1
#define HAVE_LAPACK 1

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif /* defined(PI) */

// "Vectorizes" a function call
#define LTFAT_APPLYFN(type,fn,...) do{ \
   const type list[] = {(const type)0,__VA_ARGS__}; \
   size_t len = sizeof(list)/sizeof(*list) - 1; \
   for(size_t ii=0;ii<len;ii++) \
      fn((const type)list[ii+1]); \
}while(0)

// Vectorized free
#define LTFAT_SAFEFREEALL(...) LTFAT_APPLYFN(void*,ltfat_safefree,__VA_ARGS__)

// To help muting the unused variable compiler warning
// Only works for GCC and Clang
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

#define LTFAT_MAKENAME(name,type,comp) name ## _ ## comp ## type
#define LTFAT_NAME_DOUBLE(name) LTFAT_MAKENAME(name,d,)
#define LTFAT_NAME_SINGLE(name) LTFAT_MAKENAME(name,s,)
#define LTFAT_NAME_COMPLEXDOUBLE(name) LTFAT_MAKENAME(name,d,c)
#define LTFAT_NAME_COMPLEXSINGLE(name) LTFAT_MAKENAME(name,s,c)


/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */

#ifdef MATLABFORTRAN
#define F77_FUNC(name,NAME) NAME
#else
#define F77_FUNC(name,NAME) name ## _
#endif


/* Handle Windows DLL files */
/* defined by Makefile when compiling LTFAT */
#if defined(DLL_EXPORT_SYMBOLS) && ((defined(_WIN32) || defined(__WIN32__)))
#  define LTFAT_EXTERN extern __declspec(dllexport)
#  if defined(LTFAT_DOUBLE)
#     define LTFAT_EXTERN_TOO LTFAT_EXTERN
#  else
#     define LTFAT_EXTERN_TOO
#  endif
#else
#  define LTFAT_EXTERN extern
#  define LTFAT_EXTERN_TOO LTFAT_EXTERN
#endif



/* On WinXP, gcc defines __WIN32__ */
/* On Linux, gcc defines __linux__ */

#endif /* CONFIG_H */
