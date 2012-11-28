#include <stdlib.h>
#include "ltfat.h"
#include "ltfat_time.h"

/* This timer will test only the final, unavoidable fft in the
 * computation of the rectangular DGT, to provide a baseline for
 * optimization for all the other fft routines. */


#define LTFAT_REAL double
#define LTFAT_COMPLEX ltfat_complex
#define LTFAT_NAME(name) name
#define LTFAT_FFTW(name) fftw_ ## name  


LTFAT_EXTERN LTFAT_NAME(dgt_long_plan)
LTFAT_NAME(dgt_baselinefft_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
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
   LTFAT_NAME(wfac)(g, L, 1, a, M, plan.gf);
   
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
LTFAT_NAME(dgt_baselinefft_execute)(const LTFAT_NAME(dgt_long_plan) plan)
{

   /* FFT to modulate the coefficients. */
   LTFAT_FFTW(execute)(plan.p_veryend);   

}


LTFAT_EXTERN void
LTFAT_NAME(dgt_baselinefft_done)(LTFAT_NAME(dgt_long_plan) plan)
{

   LTFAT_FFTW(destroy_plan)(plan.p_veryend);
   LTFAT_FFTW(destroy_plan)(plan.p_before);
   LTFAT_FFTW(destroy_plan)(plan.p_after);

   ltfat_free(plan.sbuf);
   ltfat_free(plan.gf);
   ltfat_free(plan.ff);
   ltfat_free(plan.cf);

}


int main( int argc, char *argv[] )
{
  ltfat_complex *f, *g, *c;
  int a, M, L, W, N, nrep, ii;
  double s0, s1;

  if (argc<5)
  {
     printf("Correct parameters: a, M, L, W, nrep\n");     
     return(1);
  }
  a = atoi(argv[1]);
  M = atoi(argv[2]);
  L = atoi(argv[3]);
  W = atoi(argv[4]);
  nrep = atoi(argv[5]);

  N=L/a;

  f  = ltfat_malloc(L*W*sizeof(ltfat_complex));
  g  = ltfat_malloc(L*W*sizeof(ltfat_complex));
  c  = ltfat_malloc(M*N*W*sizeof(ltfat_complex));
  
  dgt_long_plan plan = dgt_baselinefft_init((const ltfat_complex*)f, (const ltfat_complex*)g, L, W, a, M, c, FFTW_PATIENT);

  s0 = ltfat_time();
  for (ii=0;ii<nrep;ii++)
  {

    dgt_baselinefft_execute(plan);
    
  }
  s1 = ltfat_time();

  dgt_baselinefft_done(plan);

  printf("%i %i %i %i %f\n",a,M,L,W,(s1-s0)/nrep); 

  ltfat_free(f);
  ltfat_free(g);
  ltfat_free(c);

  return(0);
}
