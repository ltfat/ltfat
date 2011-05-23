#include "config.h"
#include "fftw3.h"
#include "ltfat.h"


LTFAT_EXTERN LTFAT_NAME(dgt_long_plan)
LTFAT_NAME(dgt_long_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			  const int L, const int W,
			  const int a, const int M, LTFAT_COMPLEX *cout,
			  unsigned flags)		      
{

   LTFAT_NAME(dgt_long_plan) plan;
   int h_m;

   plan.a=a;
   plan.M=M;
   plan.L=L;
   plan.W=W;
   const int N=L/a;
   const int b=L/M;

   plan.c=gcd(a, M,&plan.h_a, &h_m);
   const int p=a/plan.c;
   const int q=M/plan.c;
   const int d=b/p;
   plan.h_a=-plan.h_a;

   plan.sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));
   plan.cout = cout;
   plan.f    = f;

   plan.gf   = (LTFAT_COMPLEX*)ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

   plan.ff = ltfat_malloc(2*d*p*q*W*sizeof(LTFAT_REAL));
   plan.cf = ltfat_malloc(2*d*q*q*W*sizeof(LTFAT_REAL));

   /* Get factorization of window */
   LTFAT_NAME(wfac)(g, L, a, M, plan.gf);
   
  /* Create plans. In-place. */
   plan.p_veryend = 
      LTFAT_FFTW(plan_many_dft)(1, &M, N*W,
				plan.cout, NULL,
				1, M,
				plan.cout, NULL,
				1, M,
				FFTW_FORWARD, flags);
   
   plan.p_before = 
      LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)plan.sbuf,
			      (LTFAT_COMPLEX*)plan.sbuf,
			      FFTW_FORWARD, flags);
   
   plan.p_after  = 
      LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)plan.sbuf,
			      (LTFAT_COMPLEX*)plan.sbuf,
			      FFTW_BACKWARD, flags);         
   
   return plan;
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_long_execute)(const LTFAT_NAME(dgt_long_plan) plan)
{

   LTFAT_NAME(dgt_walnut_plan)(plan);
      
   /* FFT to modulate the coefficients. */
   LTFAT_FFTW(execute)(plan.p_veryend);   

}


LTFAT_EXTERN void
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan) plan)
{

   LTFAT_FFTW(destroy_plan)(plan.p_veryend);
   LTFAT_FFTW(destroy_plan)(plan.p_before);
   LTFAT_FFTW(destroy_plan)(plan.p_after);

   ltfat_free(plan.sbuf);
   ltfat_free(plan.gf);
   ltfat_free(plan.ff);
   ltfat_free(plan.cf);

}





/* -------------- Real valued signal ------------------------ */

LTFAT_EXTERN
void LTFAT_NAME(dgt_fac_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
			   const int L,
			   const int W, const int a,
			   const int M, LTFAT_COMPLEX *cout)
{

   const int N=L/a;
   
   /* Create plan. In-place. */
   LTFAT_FFTW(plan) p_veryend = 
     LTFAT_FFTW(plan_many_dft)(1, &M, N*W,
			       cout, NULL,
			       1, M,
			       cout, NULL,
			       1, M,
			       FFTW_FORWARD, FFTW_ESTIMATE);
   
   LTFAT_NAME(dgt_walnut_r)(f,gf,L,W,a,M,cout);
      
   /* FFT to modulate the coefficients. */
   LTFAT_FFTW(execute)(p_veryend);   

   LTFAT_FFTW(destroy_plan)(p_veryend);   
   
}
