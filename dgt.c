#include "config.h"
#include "fftw3.h"
#include "ltfat.h"

LTFAT_EXTERN
void LTFAT_NAME(dgt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			 const int L, const int W, const int a,
			 const int M, LTFAT_COMPLEX *cout)
{
  
  LTFAT_NAME(ltfat_plan) plan =
    LTFAT_NAME(plan_dgt_long)(f, g, L, W, a, M, cout, FFTW_ESTIMATE);
  
  LTFAT_NAME(ltfat_execute_plan)(plan);

  LTFAT_NAME(ltfat_destroy_plan)(plan);
  
}

LTFAT_EXTERN
void LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			 const int L, const int W, const int a,
			 const int M, LTFAT_COMPLEX *cout)
{
   
   const int wfs = wfacreal_size(L,a,M);

   LTFAT_COMPLEX *gf = ltfat_malloc(wfs*sizeof(LTFAT_COMPLEX));

   LTFAT_NAME(wfacreal)(g, L, a, M, gf);
   
   LTFAT_NAME(dgtreal_fac)(f, (const LTFAT_COMPLEX *)gf, L, W, a, M, cout);
  
   ltfat_free(gf);
}

