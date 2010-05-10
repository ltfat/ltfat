#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"

void LTFAT_NAME(iwfac)(const LTFAT_COMPLEX *gf, const int L,
		       const int a, const int M, LTFAT_COMPLEX *g)
{
   
   int b, N, c, d, p, q, h_a, h_m, ld3;
   
   int l, k, r, s, w;
   int rem, negrem;
   
   LTFAT_REAL scaling, *sbuf, *gfp;
   
   LTFAT_FFTW(plan) p_before;
   
   int ldf;

   const int R = 1;
   
   b=L/M;
   N=L/a;
   
   c=gcd(a, M,&h_a, &h_m);
   p=a/c;
   q=M/c;
   d=b/p;

   /* division by d is because of the way FFTW normalizes the transform. */
   scaling=1.0/sqrt(M)/d;
   
   /* for testing, set ldf=L. */
   ldf=L;

   sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));

   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
			       FFTW_BACKWARD, FFTW_MEASURE);
  
  

   ld3=c*p*q*R;
   gfp=(LTFAT_REAL*)gf;

   for (r=0;r<c;r++)
   {	
      for (w=0;w<R;w++)
      {
	 for (l=0;l<q;l++)
	 {
	    for (k=0;k<p;k++)
	    {
	       negrem=positiverem(k*M-l*a,L);
	       for (s=0;s<2*d;s+=2)
	       {	    
		  sbuf[s]   = gfp[s*ld3]*scaling;
		  sbuf[s+1] = gfp[s*ld3+1]*scaling;
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (s=0;s<d;s++)
	       {	    
		  rem = (negrem+s*p*M)%L;
		  g[r+rem+L*w][0] = sbuf[2*s];
		  g[r+rem+L*w][1] = sbuf[2*s+1];
	       }
	       gfp+=2;
	    }
	 }
      }
   }           
   
   /* Clear the work-array. */
   ltfat_free(sbuf);
}


/* This routine only returns the real valued part of the output. Use
 * this ONLY if you know beforehand that the output is real.
 */
void LTFAT_NAME(iwfac_r)(const LTFAT_COMPLEX *gf, const int L,
			 const int a, const int M, LTFAT_REAL *g)
{
   
   int b, N, c, d, p, q, h_a, h_m, ld2, ld3;
   
   int l, k, r, s, w;
   div_t domod;
   
   LTFAT_REAL scaling;
   
   LTFAT_COMPLEX *tmp_gf;
   LTFAT_FFTW(plan) p_before;
   
   const int R = 1;

   int ldf;
   
   b=L/M;
   N=L/a;
   
   c=gcd(a, M,&h_a, &h_m);
   p=a/c;
   q=M/c;
   d=b/p;

   /* division by d is because of the way FFTW normalizes the transform. */
   scaling=1.0/sqrt(M)/d;
   
   /* for testing, set ldf=L. */
   ldf=L;

   /* Create plans. In-place.
    * We must copy data to a temporary variable in order not
    * to destroy the input.
    */
   tmp_gf = (LTFAT_COMPLEX*)ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));

   /* This plan is not in-place. It copies to the work
      array. Therefore it is ok to cast away the constness of gf.*/
   p_before = LTFAT_FFTW(plan_many_dft)(1, &d, c*p*q*R,
					(LTFAT_COMPLEX *)gf, NULL,
					c*p*q*R, 1,
					tmp_gf, NULL,
					c*p*q*R, 1,
					FFTW_BACKWARD, FFTW_ESTIMATE);
     
  
   LTFAT_FFTW(execute)(p_before);	  

   ld2=p*q*R;
   ld3=c*p*q*R;

   for (w=0;w<R;w++)
   {
      for (s=0;s<d;s++)
      {
	 for (l=0;l<q;l++)
	 {
	    for (k=0;k<p;k++)
	    {
	       /*gf(k+1,l+1+q*w,r+1,s+1)=g(r+mod(k*M-l*a+s*p*M,L)+1,w+1);*/
	       /*Add L to make sure it is positive */
	       domod= div(k*M-l*a+s*p*M+L, L);
	       for (r=0;r<c;r++)
	       {	
		 /* Only copy the real part */
		 g[r+domod.rem+L*w] = tmp_gf[k+(l+q*w)*p+r*ld2+s*ld3][0]*scaling;
	       }
	    }
	 }
      }
   }           
   
   /* Clear the work-array. */
   ltfat_free(tmp_gf);
}

