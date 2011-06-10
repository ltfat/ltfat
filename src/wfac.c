#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(wfac)(const LTFAT_COMPLEX *g, const int L,
		      const int a, const int M,
		      LTFAT_COMPLEX *gf)
{
  
   int h_a, h_m, s;

   LTFAT_REAL *sbuf, *gfp;

   int rem, negrem;

   LTFAT_FFTW(plan) p_before;

   /* This is preserved for potential future use. */
   const int R = 1;
   
   const int b = L/M;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   const double sqrtM=sqrt(M);
   
   sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));
     
   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
			       FFTW_FORWARD, FFTW_MEASURE);

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
	       negrem = positiverem(k*M-l*a,L);
	       for (s=0;s<d;s++)
	       {		  
		  rem = (negrem+s*p*M)%L;
		  sbuf[2*s]   = sqrtM*g[r+rem+L*w][0];
		  sbuf[2*s+1] = sqrtM*g[r+rem+L*w][1];
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (s=0;s<2*d;s+=2)
	       {		  
		  gfp[s*ld3]  = sbuf[s];
		  gfp[s*ld3+1]= sbuf[s+1];
	       }
	       gfp+=2;
	    }
	 }
      }
   }
   
   ltfat_free(sbuf);
}


/* wfac for real valued input. */
LTFAT_EXTERN void
LTFAT_NAME(wfac_r)(const LTFAT_REAL *g, const int L,
			const int a, const int M,
			LTFAT_COMPLEX *gf)
{
  
   int h_a, h_m;   

   LTFAT_REAL *sbuf, *gfp;

   int s;
   int rem, negrem;

   const int R = 1;

   LTFAT_FFTW(plan) p_before;   
   
   const int b=L/M;
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   const double sqrtM=sqrt(M);
   
   sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));
     
   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf,
				      (LTFAT_COMPLEX*)sbuf,
				      FFTW_FORWARD, FFTW_MEASURE);

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
	       negrem = positiverem(k*M-l*a,L);
	       for (s=0;s<d;s++)
	       {		  
		  rem = (negrem+s*p*M)%L;
		  sbuf[2*s]   = sqrtM*g[r+rem+L*w];
		  sbuf[2*s+1] = 0.0;
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (s=0;s<2*d;s+=2)
	       {		  
		  gfp[s*ld3]  = sbuf[s];
		  gfp[s*ld3+1]= sbuf[s+1];
	       }
	       gfp+=2;
	    }
	 }
      }
   }
   
   ltfat_free(sbuf);
}

/* wfac for real valued input. Produces only half the output coefficients of wfac_r */
LTFAT_EXTERN void
LTFAT_NAME(wfacreal)(const LTFAT_REAL *g, const int L,
			const int a, const int M,
			LTFAT_COMPLEX *gf)
{
  
   int h_a, h_m;   

   LTFAT_REAL *gfp;

   int s;
   int rem, negrem;

   const int R = 1;

   LTFAT_FFTW(plan) p_before;   
   
   const int b=L/M;
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   /* This is a floor operation. */
   const int d2= d/2+1;

   const double sqrtM=sqrt(M);
   
   LTFAT_REAL *sbuf = (LTFAT_REAL*)ltfat_malloc(d*sizeof(LTFAT_REAL));
   LTFAT_COMPLEX *cbuf = (LTFAT_COMPLEX*)ltfat_malloc(d2*sizeof(LTFAT_COMPLEX));
     
   /* Create plan. In-place. */
   p_before = LTFAT_FFTW(plan_dft_r2c_1d)(d, sbuf, cbuf, FFTW_MEASURE);

   const int ld3=2*c*p*q*R;
   gfp=(LTFAT_REAL*)gf;
   for (int r=0;r<c;r++)
   {      
      for (int w=0;w<R;w++)
      {
	 for (int l=0;l<q;l++)
	 {
	    for (int k=0;k<p;k++)
	    {
	       negrem = positiverem(k*M-l*a,L);
	       for (s=0;s<d;s++)
	       {		  
		  rem = (negrem+s*p*M)%L;
		  sbuf[s]   = sqrtM*g[r+rem+L*w];
	       }

	       LTFAT_FFTW(execute)(p_before);	  

	       for (s=0;s<d2;s++)
	       {		  
		  gfp[s*ld3]  = cbuf[s][0];
		  gfp[s*ld3+1]= cbuf[s][1];
	       }
	       gfp+=2;
	    }
	 }
      }
   }
   
   ltfat_free(sbuf);
}
