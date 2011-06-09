#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(ufilterbank_fft)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			    const int L, const int gl,
			    const int W, const int a, const int M,
			    LTFAT_COMPLEX *cout)
{

   /* ----- Initialization ------------ */

   const int N=L/a;

   LTFAT_COMPLEX *gwork = ltfat_malloc(L*M*sizeof(LTFAT_COMPLEX));

   LTFAT_COMPLEX *work = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

   LTFAT_FFTW(plan) plan_g = 
      LTFAT_FFTW(plan_many_dft)(1, &L, M,
				gwork, NULL,
				1, L,
				gwork, NULL,
				1, L,
				FFTW_FORWARD, FFTW_ESTIMATE);

      LTFAT_FFTW(plan_dft_1d)(L, gwork, gwork,
			      FFTW_FORWARD, FFTW_ESTIMATE);
   
   LTFAT_FFTW(plan) plan_w =
       LTFAT_FFTW(plan_dft_1d)(L, work, work,
			      FFTW_FORWARD, FFTW_ESTIMATE);

   LTFAT_FFTW(plan) plan_c = 
      LTFAT_FFTW(plan_many_dft)(1, &N, M*W,
				cout, NULL,
				1, N,
				cout, NULL,
				1, N,
				FFTW_BACKWARD, FFTW_ESTIMATE);

   const LTFAT_REAL scalconst = 1.0/L;

   /* ----- Main -------------------------- */

   /* Extend g and copy to work buffer */
   for (int m=0; m<M; m++)
   {
      LTFAT_NAME(fir2long_c)(g+m*gl, gl, L, gwork+m*L);
   }

   LTFAT_FFTW(execute)(plan_g);

   for (int w=0; w<W; w++)
   {      
      memcpy(work,f+L*w,sizeof(LTFAT_COMPLEX)*L);
      LTFAT_FFTW(execute)(plan_w);
	 
      for (int m=0; m<M; m++)
      {	 
	 for (int n=0; n<N; n++)
	 {
	    cout[n+m*N+w*N*M][0]=0.0;
	    cout[n+m*N+w*N*M][1]=0.0;
	    for (int k=0; k<a; k++)
	    {
	       const int l=n+k*N;
	       const LTFAT_REAL tmp0 = work[l][0]*gwork[l+m*L][0]-work[l][1]*gwork[l+m*L][1];
	       const LTFAT_REAL tmp1 = work[l][0]*gwork[l+m*L][1]+work[l][1]*gwork[l+m*L][0];
	       cout[n+m*N+w*N*M][0]+=tmp0*scalconst;
	       cout[n+m*N+w*N*M][1]+=tmp1*scalconst;
	    }
	 }
      }
   }

   
   LTFAT_FFTW(execute)(plan_c);

   

   ltfat_free(work);
   ltfat_free(gwork);

}

