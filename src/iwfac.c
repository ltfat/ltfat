#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(iwfac)(const LTFAT_COMPLEX *gf, const int L,
		       const int a, const int M, LTFAT_COMPLEX *g)
{
   
   int h_a, h_m;
   
   int rem, negrem;
   
   LTFAT_REAL scaling, *sbuf, *gfp;
   
   LTFAT_FFTW(plan) p_before;
   
   const int R = 1;
   
   const int b=L/M;
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   /* division by d is because of the way FFTW normalizes the transform. */
   scaling=1.0/sqrt(M)/d;
   
   sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));

   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
			       FFTW_BACKWARD, FFTW_MEASURE);
  
  

   const int ld3=c*p*q*R;
   gfp=(LTFAT_REAL*)gf;

   for (int r=0;r<c;r++)
   {	
      for (int w=0;w<R;w++)
      {
	 for (int l=0;l<q;l++)
	 {
	    for (int k=0;k<p;k++)
	    {
	       negrem=positiverem(k*M-l*a,L);
	       for (int s=0;s<2*d;s+=2)
	       {	    
		  sbuf[s]   = gfp[s*ld3]*scaling;
		  sbuf[s+1] = gfp[s*ld3+1]*scaling;
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (int s=0;s<d;s++)
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



LTFAT_EXTERN void
LTFAT_NAME(iwfac_r)(const LTFAT_COMPLEX *gf, const int L,
		       const int a, const int M, LTFAT_REAL *g)
{
   
   int h_a, h_m;
   
   int rem, negrem;
   
   LTFAT_REAL scaling, *sbuf, *gfp;
   
   LTFAT_FFTW(plan) p_before;
   
   const int R = 1;
   
   const int b=L/M;
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   /* division by d is because of the way FFTW normalizes the transform. */
   scaling=1.0/sqrt(M)/d;
   
   sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));

   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
			       FFTW_BACKWARD, FFTW_MEASURE);
  
  

   const int ld3=c*p*q*R;
   gfp=(LTFAT_REAL*)gf;

   for (int r=0;r<c;r++)
   {	
      for (int w=0;w<R;w++)
      {
	 for (int l=0;l<q;l++)
	 {
	    for (int k=0;k<p;k++)
	    {
	       negrem=positiverem(k*M-l*a,L);
	       for (int s=0;s<2*d;s+=2)
	       {	    
		  sbuf[s]   = gfp[s*ld3]*scaling;
		  sbuf[s+1] = gfp[s*ld3+1]*scaling;
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (int s=0;s<d;s++)
	       {	    
		  rem = (negrem+s*p*M)%L;
		  g[r+rem+L*w] = sbuf[2*s];
	       }
	       gfp+=2;
	    }
	 }
      }
   }           
   
   /* Clear the work-array. */
   ltfat_free(sbuf);
}


LTFAT_EXTERN void
LTFAT_NAME(iwfacreal)(const LTFAT_COMPLEX *gf, const int L,
		       const int a, const int M, LTFAT_REAL *g)
{
   
   int h_a, h_m;
      
   LTFAT_FFTW(plan) p_before;   
   
   const int R = 1;
   
   const int b=L/M;
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   /* This is a floor operation. */
   const int d2= d/2+1;

   /* division by d is because of the way FFTW normalizes the transform. */
   const LTFAT_REAL scaling=1.0/sqrt(M)/d;
   
   LTFAT_REAL    *sbuf = ltfat_malloc( d*sizeof(LTFAT_REAL));
   LTFAT_COMPLEX *cbuf = ltfat_malloc(d2*sizeof(LTFAT_COMPLEX));

   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_c2r_1d)(d, cbuf, sbuf, FFTW_MEASURE);

   const int ld3=c*p*q*R;

   /* Advancing pointer: Runs through array pointing out the base for the strided operations. */
   const LTFAT_COMPLEX *gfp = gf;

   for (int r=0;r<c;r++)
   {	
      for (int w=0;w<R;w++)
      {
	 for (int l=0;l<q;l++)
	 {
	    for (int k=0;k<p;k++)
	    {
	       const int negrem=positiverem(k*M-l*a,L);
	       for (int s=0;s<d2;s++)
	       {	    
		  cbuf[s][0] = gfp[s*ld3][0]*scaling;
		  cbuf[s][1] = gfp[s*ld3][1]*scaling;
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (int s=0;s<d;s++)
	       {	    
		  g[r+(negrem+s*p*M)%L+L*w] = sbuf[s];
	       }
	       gfp++;
	    }
	 }
      }
   }           
   
   /* Clear the work-array. */
   ltfat_free(cbuf);
   ltfat_free(sbuf);
}



