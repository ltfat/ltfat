#ifndef LTFAT_H
#define LTFAT_H 1

#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
//#include "config.h"
#include "cblas.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

/*
SAFE SIGNED AND UNSIGNED INTEGERS
Currently not used anywhere.
*/
typedef size_t    Lsize;
typedef size_t    Lindex;
typedef ptrdiff_t Lint;

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



  /* -------- Define routines that do not change between single/double-- */

int gcd(const int a, const int b, int *r, int *s );


void* ltfat_malloc (size_t n);


void* ltfat_calloc (size_t nmemb, size_t size);


void* ltfat_realloc (void *ptr, size_t n);


void* ltfat_realloc_and_copy (void *ptr, size_t nold, size_t nnew);


void  ltfat_free(void *ptr);
void  ltfat_safefree(void *ptr);

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

/*
MACRO-FU
*/

// "Vectorizes" function call
#define LTFAT_APPLYFN(type,fn,...) do{ \
   type* list = (type[]) {(type)0,__VA_ARGS__}; \
   size_t len = sizeof(list)/sizeof(*list) - 1; \
   for(size_t ii=0;ii<len;ii++) \
      fn((type)list[ii+1]); \
}while(0)

#define LTFAT_SAFEFREEALL(...) LTFAT_APPLYFN(void*,ltfat_safefree,__VA_ARGS__)

/* END_C_DECLS */

#endif /* !LTFAT_H */
