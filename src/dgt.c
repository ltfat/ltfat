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

