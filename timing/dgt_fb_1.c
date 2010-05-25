#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "dgt.h"
#include "tfutil.h"
#include "ltfat_time.h"


static inline int complexprod(ltfat_complex *c, ltfat_complex a,ltfat_complex b)
{
#ifdef HAVE_COMPLEX_H
  (*c)=a*b;
#else
  
  (*c)[0]=a[0]*b[0]-a[1]*b[1];
  (*c)[1]=a[1]*b[0]+a[0]*b[1];

#endif
  return (0);
}

/*  Basic version of the filter bank algorithm. The modulation is
    applied in a separate step in the end. 

    This algorithm DOES NOT HANDLE MULTIWINDOWS CORRECTLY
*/

void dgt_fb(ltfat_complex *f, ltfat_complex *g,
	    const int L, const int gl,
	    const int W, const int R, const int a, const int M, 
	    ltfat_complex *cout, int dotime)
{
   /*  --------- initial declarations -------------- */

   int b, N;
   int glh, glh_d_a;

   int k, l, m, n, r, w, rem;

   ltfat_complex *gw, *fw;

   fftw_plan p_veryend;

   ltfat_complex *coefsum;

   double st0, st1, st2, st3;
   
   /*  ----------- calculation of parameters and plans -------- */
   
   if (dotime)
   {
      st0=ltfat_time();
   }

   b=L/M;
   N=L/a;

   /* This is a floor operation. */
   glh=gl/2;

   /* This is a ceil operation. */
   glh_d_a=(int)ceil((glh*1.0)/(a));
   
   gw = (ltfat_complex*)ltfat_malloc(gl*R*sizeof(ltfat_complex));
   fw = (ltfat_complex*)ltfat_malloc(gl*W*R*sizeof(ltfat_complex));

   /* Create plan. In-place. */
   p_veryend = fftw_plan_many_dft(1, &M, N*R*W,
				  cout, NULL,
				  1, M,
				  cout, NULL,
				  1, M,
				  FFTW_FORWARD, FFTW_OPTITYPE);
      
   if (dotime)
   {
      st1=ltfat_time();
      printf("DGT_FB_1: Planning phase %f\n",st1-st0);
   }

   /*  ---------- main body ----------- */
  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */
   
   for (r=0;r<R;r++)
   {
      for (l=0;l<glh;l++)
      {
#ifdef HAVE_COMPLEX_H
	 gw[l+gl*r]=conj(g[l+(gl-glh)+gl*r]);
#else
	 gw[l+gl*r][0]=g[l+(gl-glh)+gl*r][0];
	 gw[l+gl*r][1]=-g[l+(gl-glh)+gl*r][1];
#endif
      }
      for (l=glh;l<gl;l++)
      {
#ifdef HAVE_COMPLEX_H
	 gw[l+gl*r]=conj(g[l-glh+gl*r]);
#else
	 gw[l+gl*r][0]=g[l-glh+gl*r][0];
	 gw[l+gl*r][1]=-g[l-glh+gl*r][1];
#endif
      }
   }

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (n=0; n<glh_d_a; n++)
   {
      
      /* Matlab code for this:
       * fpart=[f(L-(glh-n*a)+1:L,:);...
       * f(1:gl-(glh-n*a),:)];
       * 
       * fg=fpart.*gw;
       */
      for (r=0;r<R;r++)
      {
	 for (w=0;w<W;w++)
	 {
	    for (l=0;l<glh-n*a;l++)
	    {
	      complexprod(fw+l+gl*w+gl*W*r,f[L-(glh-n*a)+l+L*w],gw[l+r*gl]);
	    }
	    for (l=glh-n*a;l<gl;l++)
	    {
	      complexprod(fw+l+gl*w+gl*W*r,f[l-(glh-n*a)+L*w],gw[l+r*gl]);
	    }
	 }	 
      }

      /* Do the sum (decimation in frequency, Poisson summation) */      
      /*      
       *coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
       */
      for (m=0;m<M;m++)
      {
	 /* Place coefficient correctly such that we get a
	  * frequency invariant system. 
	  */
	 rem = positiverem(m+(n*a-glh), M);
	 coefsum=cout+n*M+rem;
	 for (w=0;w<W;w++)
	 {	 	
#ifdef HAVE_COMPLEX_H  
	    coefsum[w*M*N]=0.0;
#else
	    coefsum[w*M*N][0]=0.0;
	    coefsum[w*M*N][1]=0.0;
#endif
	    for (k=0;k<gl/M;k++)
	    {
#ifdef HAVE_COMPLEX_H
	       coefsum[w*M*N]+= fw[m+k*M+gl*w];
#else
	       coefsum[w*M*N][0]+= fw[m+k*M+gl*w][0];
	       coefsum[w*M*N][1]+= fw[m+k*M+gl*w][1];
#endif
	    }
	 }
      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {

      /*
	fg=f(n*a-glh+1:n*a-glh+gl,:).*gw;
	
	% Do the sum (decimation in frequency, Poisson summation)
	coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
	end;
      */
      for (r=0;r<R;r++)
      {
	 for (w=0;w<W;w++)
	 {
	    for (l=0;l<gl;l++)
	    {
	      complexprod(fw+l+gl*w+gl*W*r,f[l+n*a-glh+L*w],gw[l+r*gl]);
	    }
	 }	 
      }

      /* Do the sum (decimation in frequency, Poisson summation) */      
      /*      
       *coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
       */
      for (m=0;m<M;m++)
      {
	 /* Place coefficient correctly such that we get a
	  * frequency invariant system.
	  */
	 rem = positiverem(m+(n*a-glh), M);
	 coefsum=cout+n*M+rem;
	 for (w=0;w<W;w++)
	 {	 	
#ifdef HAVE_COMPLEX_H  
	    coefsum[w*M*N]=0.0;
#else
	    coefsum[w*M*N][0]=0.0;
	    coefsum[w*M*N][1]=0.0;
#endif
	    for (k=0;k<gl/M;k++)
	    {
#ifdef HAVE_COMPLEX_H
	       coefsum[w*M*N]+= fw[m+k*M+gl*w];
#else
	       coefsum[w*M*N][0]+= fw[m+k*M+gl*w][0];
	       coefsum[w*M*N][1]+= fw[m+k*M+gl*w][1];
#endif
	    }
	 }
      }
   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      /*
	fpart=[f((n*a-glh)+1:L,:);... %   L-n*a+glh) elements
	f(1:n*a-glh+gl-L,:)];  %  gl-L+n*a-glh elements
	fg=fpart.*gw;
	
	% Do the sum (decimation in frequency, Poisson summation)
	coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);      
      */
      for (r=0;r<R;r++)
      {
	 for (w=0;w<W;w++)
	 {
	    for (l=0;l<L-n*a+glh;l++)
	    {
	      complexprod(fw+l+gl*w+gl*W*r,f[n*a-glh+l+L*w],gw[l+r*gl]);
	    }
	    for (l=L-n*a+glh;l<gl;l++)
	    {
	      complexprod(fw+l+gl*w+gl*W*r,f[l-(L-n*a+glh)+L*w],gw[l+r*gl]);
	    }
	 }	 
      }


      /* Do the sum (decimation in frequency, Poisson summation) */      
      /*      
       *coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
       */
      for (m=0;m<M;m++)
      {
	 /* Place coefficient correctly such that we get a
	  * frequency invariant system. Add 2*L to make sure it is always
	  * positive.
	  */
	 rem = positiverem(m+(n*a-glh), M);
	 coefsum=cout+n*M+rem;
	 for (w=0;w<W;w++)
	 {	 	
#ifdef HAVE_COMPLEX_H
	    coefsum[w*M*N]=0.0;
#else
	    coefsum[w*M*N][0]=0.0;
	    coefsum[w*M*N][1]=0.0;
#endif
	    for (k=0;k<gl/M;k++)
	    {
#ifdef HAVE_COMPLEX_H
	       coefsum[w*M*N]+= fw[m+k*M+gl*w];
#else
	       coefsum[w*M*N][0]+= fw[m+k*M+gl*w][0];
	       coefsum[w*M*N][1]+= fw[m+k*M+gl*w][1];
#endif
	    }
	 }
      }

   }

   if (dotime)
   {
      st2=ltfat_time();
      printf("DGT_FB_1: Window application and poisson summing done %f\n",st2-st1);
   }


   /* FFT to modulate the coefficients. */
   fftw_execute(p_veryend);   

   if (dotime)
   {
      st3=ltfat_time();
      printf("DGT_FB_1: Final FFT  %f\n",st3-st2);
      printf("DGT_FB_1: Total time %f\n",st3-st0);   
   }

    /* -----------  Clean up ----------------- */   
   ltfat_free(gw);
   ltfat_free(fw);    

   fftw_destroy_plan(p_veryend);
}
