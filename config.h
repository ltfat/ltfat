#ifndef CONFIG_H
#define CONFIG_H 1


#ifndef __cplusplus
#include <complex.h>
#endif //__cplusplus

#define HAVE_BLAS 1
#define HAVE_LAPACK 1

#include "fftw3.h"


static inline int ltfat_round(double x)
{
    if (x < 0.0)
        return (int)(x - 0.5);
    else
        return (int)(x + 0.5);
}

static inline int positiverem(int a,int b)
{
    const int c = a%b;
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
