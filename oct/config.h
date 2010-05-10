/* This file should contain the configuration parameters
   necessary for the Oct-compilation. */

#ifndef CONFIG_H
#define CONFIG_H 1

#include "fftw3.h"
#include "ltfat.h"

#define FFTW_OPTITYPE FFTW_ESTIMATE

typedef fftw_complex  ltfat_complex;
typedef fftwf_complex ltfat_scomplex;

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

#endif /* CONFIG_H */
