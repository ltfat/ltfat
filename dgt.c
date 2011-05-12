#include <string.h>
#include "config.h"
#include "fftw3.h"
#include "ltfat.h"

LTFAT_EXTERN
void LTFAT_NAME(dgt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			 const int L, const int W, const int a,
			 const int M, LTFAT_COMPLEX *cout)
{
  
  LTFAT_NAME(dgt_long_plan) plan =
    LTFAT_NAME(dgt_long_init)(f, g, L, W, a, M, cout, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgt_long_execute)(plan);

  LTFAT_NAME(dgt_long_done)(plan);
  
}

LTFAT_EXTERN
void LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			 const int L, const int W, const int a,
			 const int M, LTFAT_COMPLEX *cout)
{
 
  LTFAT_NAME(dgtreal_long_plan) plan =
    LTFAT_NAME(dgtreal_long_init)(f, g, L, W, a, M, cout, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgtreal_long_execute)(plan);

  LTFAT_NAME(dgtreal_long_done)(plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
		     const int L, const int gl,
		     const int W,  const int a, const int M, 
		     LTFAT_COMPLEX *cout)
{
 
  LTFAT_NAME(dgt_fb_plan) plan =
    LTFAT_NAME(dgt_fb_init)(g, gl, a, M, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgt_fb_execute)(plan, f, L, W, cout);

  LTFAT_NAME(dgt_fb_done)(plan);

}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
	      const int L, const int gl,
	      const int W, const int a, const int M, 
	      LTFAT_COMPLEX *cout)
{
 
  LTFAT_NAME(dgtreal_fb_plan) plan =
    LTFAT_NAME(dgtreal_fb_init)(g, gl, a, M, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgtreal_fb_execute)(plan, f, L, W, cout);

  LTFAT_NAME(dgtreal_fb_done)(plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
		    const int L, const int gl,
		    const int W, const int a, const int M, const int bl, 
		    LTFAT_COMPLEX *cout)
{
  LTFAT_NAME(dgt_ola_plan) plan =
    LTFAT_NAME(dgt_ola_init)(g, gl, W, a, M, bl, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgt_ola_execute)(plan, f, L, cout);

  LTFAT_NAME(dgt_ola_done)(plan);

}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola)(const LTFAT_REAL *f, const LTFAT_REAL *g,
		    const int L, const int gl,
		    const int W, const int a, const int M, const int bl, 
		    LTFAT_COMPLEX *cout)
{
  LTFAT_NAME(dgtreal_ola_plan) plan =
    LTFAT_NAME(dgtreal_ola_init)(g, gl, W, a, M, bl, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgtreal_ola_execute)(plan, f, L, cout);

  LTFAT_NAME(dgtreal_ola_done)(plan);

}


