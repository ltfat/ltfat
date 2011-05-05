#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"


LTFAT_NAME(dgtreal_long_plan)
LTFAT_NAME(dgtreal_long_init)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			  const int L, const int W,
			  const int a, const int M, LTFAT_COMPLEX *cout,
			  unsigned flags)		      
{

   LTFAT_NAME(dgtreal_long_plan) plan;
   int h_m;

   plan.a=a;
   plan.M=M;
   plan.L=L;
   plan.W=W;
   const int N=L/a;


   plan.c=gcd(a, M,&plan.h_a, &h_m);
   const int b=L/M;
   const int p=a/plan.c;
   const int d=b/p;
   plan.h_a=-plan.h_a;

   const int M2=M/2+1;
   const int d2=d/2+1;

   plan.sbuf = ltfat_malloc( d*sizeof(LTFAT_REAL));
   plan.cbuf = ltfat_malloc(d2*sizeof(LTFAT_COMPLEX));
   plan.cout = cout;
   plan.f    = f;

   const int wfs = wfacreal_size(L,a,M);

   plan.gf   = (LTFAT_COMPLEX*)ltfat_malloc(wfs*sizeof(LTFAT_COMPLEX));

   plan.cwork = (LTFAT_REAL*)ltfat_malloc(M*N*W*sizeof(LTFAT_REAL));

   /* Get factorization of window */
   LTFAT_NAME(wfacreal)(g, L, a, M, plan.gf);
   
  /* Create plans. In-place. */
   plan.p_veryend = 
      LTFAT_FFTW(plan_many_dft_r2c)(1, &plan.M, N*W,
                                  plan.cwork, NULL,
                                  1, M,
                                  cout, NULL,
				  1, M2,
				  flags);

   plan.p_before = 
      LTFAT_FFTW(plan_dft_r2c_1d)(d, plan.sbuf, plan.cbuf, flags);
   
   plan.p_after  = 
      LTFAT_FFTW(plan_dft_c2r_1d)(d, plan.cbuf, plan.sbuf, flags);         
   
   return plan;
}




void LTFAT_NAME(dgtreal_long_execute)(const LTFAT_NAME(dgtreal_long_plan) plan)
{

   LTFAT_NAME(dgtreal_walnut_plan)(plan);
      
   /* FFT to modulate the coefficients. */
   LTFAT_FFTW(execute)(plan.p_veryend);   

}


void LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan) plan)
{

   LTFAT_FFTW(destroy_plan)(plan.p_veryend);
   LTFAT_FFTW(destroy_plan)(plan.p_before);
   LTFAT_FFTW(destroy_plan)(plan.p_after);

   ltfat_free(plan.sbuf);
   ltfat_free(plan.cbuf);
   ltfat_free(plan.cwork);
   ltfat_free(plan.gf);
}






/* LTFAT_EXTERN void */
/* LTFAT_NAME(dgtreal_fac_deprecated)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf, */
/* 			     const int L,const int W, const int a,	      */
/* 			     const int M, LTFAT_COMPLEX *cout) */
/* { */
/*    LTFAT_FFTW(plan) p_veryend; */
/*    LTFAT_REAL *cwork; */

/*    /\*  ----------- calculation of parameters and plans -------- *\/ */
      
/*    const int N=L/a; */

/*    /\* This is a floor operation. *\/ */
/*    const int M2= M/2+1; */

/*    cwork = (LTFAT_REAL*)ltfat_malloc(M*N*W*sizeof(LTFAT_REAL)); */
   
/*    /\* Create plan. In-place. *\/ */
/*    p_veryend = LTFAT_FFTW(plan_many_dft_r2c)(1, &M, N*W, */
/*                                   cwork, NULL, */
/*                                   1, M, */
/*                                   cout, NULL, */
/* 				  1, M2, */
/* 				  FFTW_ESTIMATE); */

/*    LTFAT_NAME(dgtreal_walnut)(f,gf,L,W,a,M,cwork); */
      
/*    /\* FFT to modulate the coefficients. *\/ */
/*    LTFAT_FFTW(execute)(p_veryend);    */

/*    LTFAT_FFTW(destroy_plan)(p_veryend); */

/*    ltfat_free(cwork); */
   
/* } */
