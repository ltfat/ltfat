#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"

void LTFAT_NAME(dgtreal_fac)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
			     const int L,const int W, const int a,	     
			     const int M, LTFAT_COMPLEX *cout)
{
   LTFAT_FFTW(plan) p_veryend;
   LTFAT_REAL *cwork;

   /*  ----------- calculation of parameters and plans -------- */
      
   const int N=L/a;

   /* This is a floor operation. */
   const int M2= M/2+1;

   cwork = (LTFAT_REAL*)ltfat_malloc(M*N*W*sizeof(LTFAT_REAL));
   
   /* Create plan. In-place. */
   p_veryend = LTFAT_FFTW(plan_many_dft_r2c)(1, &M, N*W,
                                  cwork, NULL,
                                  1, M,
                                  cout, NULL,
				  1, M2,
				  FFTW_ESTIMATE);

   LTFAT_NAME(dgtreal_walnut)(f,gf,L,W,a,M,cwork);
      
   /* FFT to modulate the coefficients. */
   LTFAT_FFTW(execute)(p_veryend);   

   LTFAT_FFTW(destroy_plan)(p_veryend);

   ltfat_free(cwork);
   
}
