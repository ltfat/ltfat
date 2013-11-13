#ifndef LTFAT_H
#define LTFAT_H 1

#include <stdlib.h>
#include <stddef.h>
#include "fftw3.h"
#include "cblas.h"

// SAFE SIGNED AND UNSIGNED INTEGERS
typedef size_t    ltfatSize;
typedef size_t    ltfatIndex;
typedef ptrdiff_t ltfatDiff;

/* Allow using arbitrary complex types. Assumes identical memory layout. */
#ifndef LTFAT_USER_COMPLEX
typedef fftw_complex  ltfat_complex;
typedef fftwf_complex ltfat_scomplex;
#endif

#define LTFAT_MAKENAME(name,type,comp) name ## _ ## comp ## type
#define LTFAT_NAME_DOUBLE(name) LTFAT_MAKENAME(name,d,)
#define LTFAT_NAME_SINGLE(name) LTFAT_MAKENAME(name,s,)
#define LTFAT_NAME_COMPLEXDOUBLE(name) LTFAT_MAKENAME(name,d,c)
#define LTFAT_NAME_COMPLEXSINGLE(name) LTFAT_MAKENAME(name,s,c)

/* BEGIN_C_DECLS */

#ifdef __cplusplus
extern "C"
{
#endif


/* Handle Windows DLL files, not used */

#if defined(LTFAT_DLL_NEVERUSED) && (defined(_WIN32) || defined(__WIN32__))
#  if defined(COMPILING_LTFAT) /* defined by Makefile when compiling LTFAT */
#    define LTFAT_EXTERN extern __declspec(dllexport)
#    if defined(LTFAT_DOUBLE)
#      define LTFAT_EXTERN_NOTYPE extern __declspec(dllexport)
#    else
#      define LTFAT_EXTERN_NOTYPE extern
#    endif
#  else /* user is calling LTFAT; import symbol */
#    define LTFAT_EXTERN extern __declspec(dllimport)
#  endif
#else
#  define LTFAT_EXTERN extern
#  define LTFAT_EXTERN_NOTYPE extern
#endif

  /* -------- Define the double precision routines ----- */

#define LTFAT_H_REAL double
#define LTFAT_H_COMPLEX fftw_complex
//#define LTFAT_H_COMPLEX _Complex double
#define LTFAT_H_COMPLEXH double _Complex
#define LTFAT_H_TYPE LTFAT_H_REAL
#define LTFAT_H_NAME(name) LTFAT_NAME_DOUBLE(name)
#define LTFAT_H_NAME_REAL(name) LTFAT_NAME_DOUBLE(name)
#define LTFAT_H_NAME_COMPLEX(name) LTFAT_NAME_COMPLEXDOUBLE(name)
#define LTFAT_H_FFTW(name) fftw_ ## name

#include "ltfat_typeindependent.h"
#undef LTFAT_H_COMPLEX
#define LTFAT_H_COMPLEX LTFAT_H_COMPLEXH
#include "ltfat_typecomplexindependent.h"

#undef LTFAT_H_TYPE
#undef LTFAT_H_NAME
#define LTFAT_H_TYPE LTFAT_H_COMPLEX
#define LTFAT_H_NAME(name) LTFAT_NAME_COMPLEXDOUBLE(name)

#include "ltfat_typecomplexindependent.h"

#undef LTFAT_H_REAL
#undef LTFAT_H_COMPLEX
#undef LTFAT_H_COMPLEXH
#undef LTFAT_H_NAME
#undef LTFAT_H_FFTW
#undef LTFAT_H_TYPE
#undef LTFAT_H_NAME_COMPLEX


  /* -------- Define the single precision routines ----- */

#define LTFAT_H_REAL float
#define LTFAT_H_COMPLEX fftwf_complex
//#define LTFAT_H_COMPLEX _Complex float
#define LTFAT_H_COMPLEXH float _Complex
#define LTFAT_H_NAME(name) LTFAT_NAME_SINGLE(name)
#define LTFAT_H_NAME_COMPLEX(name) LTFAT_NAME_COMPLEXSINGLE(name)
#define LTFAT_H_FFTW(name) fftwf_ ## name
#define LTFAT_H_TYPE LTFAT_H_REAL


#include "ltfat_typeindependent.h"
#undef LTFAT_H_COMPLEX
#define LTFAT_H_COMPLEX LTFAT_H_COMPLEXH

#include "ltfat_typecomplexindependent.h"

#undef LTFAT_H_TYPE
#undef LTFAT_H_NAME
#define LTFAT_H_TYPE LTFAT_H_COMPLEX
#define LTFAT_H_NAME(name) LTFAT_NAME_COMPLEXSINGLE(name)


#include "ltfat_typecomplexindependent.h"

#undef LTFAT_H_REAL
#undef LTFAT_H_COMPLEX
#undef LTFAT_H_COMPLEXH
#undef LTFAT_H_TYPE
#undef LTFAT_H_NAME
#undef LTFAT_H_FFTW
#undef LTFAT_H_NAME_COMPLEX


  /* -------- Define routines that do not change between single/double-- */
//LTFAT_EXTERN_NOTYPE
int gcd(const int a, const int b, int *r, int *s );

//LTFAT_EXTERN_NOTYPE
void* ltfat_malloc (size_t n);

//LTFAT_EXTERN_NOTYPE
void* ltfat_calloc (size_t nmemb, size_t size);

//LTFAT_EXTERN_NOTYPE
void* ltfat_realloc (void *ptr, size_t n);

//LTFAT_EXTERN_NOTYPE
void* ltfat_realloc_and_copy (void *ptr, size_t nold, size_t nnew);

//LTFAT_EXTERN_NOTYPE
void  ltfat_free(void *ptr);

//LTFAT_EXTERN_NOTYPE
void fftindex(const int N, int *indexout);

//LTFAT_EXTERN_NOTYPE
int makelarger(const int L, const int K);

//LTFAT_EXTERN_NOTYPE
int int_max(const int a, const int b);

//LTFAT_EXTERN_NOTYPE
int int_min(const int a, const int b);

//LTFAT_EXTERN_NOTYPE
int lcm(const int a, const int b);

//LTFAT_EXTERN_NOTYPE
void gabimagepars(const int Ls, const int x, const int y,
		  int *a, int *M, int *L, int *N, int *Ngood);

//LTFAT_EXTERN_NOTYPE
int wfacreal_size(const int L, const int a, const int M);


size_t nextPow2_st(size_t x);

size_t nextfastfft(size_t x);



size_t max_st(const size_t a, const size_t b);
size_t min_st(const size_t a, const size_t b);

ptrdiff_t max_pt(const ptrdiff_t a, const ptrdiff_t b);
ptrdiff_t min_pt(const ptrdiff_t a, const ptrdiff_t b);


#ifdef __cplusplus
}  // extern "C"
#endif


/* END_C_DECLS */

#endif /* !LTFAT_H */
