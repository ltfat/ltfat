#ifndef CONFIG_H
#define CONFIG_H 1

#define HAVE_BLAS 1
#define HAVE_LAPACK 1

#include <stdlib.h>
#include "stddef.h"

#include "fftw3.h"

#ifdef LTFAT_DOUBLE
#define LTFAT_COMPLEX fftw_complex
#define LTFAT_REAL double
#define LTFAT_NAME(name) name  
#define LTFAT_FFTW(name) fftw_ ## name  
#endif

#ifdef LTFAT_SINGLE
#define LTFAT_COMPLEX fftwf_complex
#define LTFAT_REAL float
#define LTFAT_NAME(name) s ## name
#define LTFAT_FFTW(name) fftwf_ ## name  
#endif


static inline int ltfat_round(double x)
{
  return (int)(x+.5); 
}

static inline int positiverem(int a,int b)
{
  int c;
  c=a%b;
  return(c<0 ? c+b : c);
}

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */

#ifdef MATLABFORTRAN
#define F77_FUNC(name,NAME) NAME 
#else
#define F77_FUNC(name,NAME) name ## _
#endif

 /* On WinXP, gcc defines __WIN32__ */
 /* On Linux, gcc defines __linux__ */

#endif /* CONFIG_H */
