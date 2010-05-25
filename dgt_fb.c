#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "ltfat.h"

/* The following macro adds the coefficients together performing the
 * last part of the Poisson summation, executes the FFT on the summed
 * coefficients, and places the coefficients in the output array.
 *
 * The first summation is done in that peculiar way to obtain the
 * correct phase for a frequency invariant Gabor transform. Summing
 * them directly would lead to a time invariant (phase-locked) Gabor
 * transform.
 *
 * The macro is called in three different places in the dgt_fb function.
 */
#define THE_SUM { \
for (m=0;m<M;m++) \
{ \
   rem = 2*positiverem(m+(n*a-glh), M); \
   sbuf[rem]=0.0; \
   sbuf[rem+1]=0.0; \
   fbd=fw+2*m; \
   for (k=0;k<gl/M;k++) \
   { \
      sbuf[rem]  += fbd[0]; \
      sbuf[rem+1]+= fbd[1]; \
      fbd+=2*M; \
   } \
} \
\
 LTFAT_FFTW(execute)(p_small);			\
\
coefsum=(LTFAT_REAL*)cout+2*(n*M+r*M*N+w*M*N*R); \
for (m=0;m<M;m++) \
{ \
   coefsum[2*m] = sbuf[2*m]; \
   coefsum[2*m+1] = sbuf[2*m+1]; \
}} 

#define THE_SUM_REAL { \
for (m=0;m<M;m++) \
{ \
   rem = positiverem(m+(n*a-glh), M); \
   sbuf_in[rem]=0.0; \
   fbd=fw+m; \
   for (k=0;k<gl/M;k++) \
   { \
     sbuf_in[rem]+=(*fbd);			\
      fbd+=M; \
   } \
} \
\
 LTFAT_FFTW(execute)(p_small);			\
\
coefsum=(LTFAT_REAL*)cout+2*(n*M2+r*M2*N+w*M2*N*R); \
for (m=0;m<M2;m++) \
{ \
   coefsum[2*m]   = sbuf_out[m][0]; \
   coefsum[2*m+1] = sbuf_out[m][1]; \
}}


void LTFAT_NAME(dgt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
	      const int L, const int gl,
	      const int W, const int a, const int M, 
	      LTFAT_COMPLEX *cout)
{
   /*  --------- initial declarations -------------- */

   int b, N;
   int glh, glh_d_a;

   int k, l, m, n, r, w, rem;

   LTFAT_COMPLEX *gw;

   LTFAT_FFTW(plan) p_small;

   LTFAT_COMPLEX *gb;
   LTFAT_REAL *sbuf,*coefsum, *fbd, *fw;
   
   const int R = 1;

   /*  ----------- calculation of parameters and plans -------- */
   
   b=L/M;
   N=L/a;

   /* This is a floor operation. */
   glh=gl/2;

   /* This is a ceil operation. */
   glh_d_a=(int)ceil((glh*1.0)/(a));
   
   gw   = (LTFAT_COMPLEX*)ltfat_malloc(gl*R*sizeof(LTFAT_COMPLEX));
   fw   = (LTFAT_REAL*)ltfat_malloc(2*gl*sizeof(LTFAT_REAL));
   sbuf = (LTFAT_REAL*)ltfat_malloc(2*M*sizeof(LTFAT_REAL));

   /* Create plan. In-place. */
   p_small = LTFAT_FFTW(plan_dft_1d)(M, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
			      FFTW_FORWARD, FFTW_MEASURE);
   
   /*  ---------- main body ----------- */
  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */
   
   for (r=0;r<R;r++)
   {
      for (l=0;l<glh;l++)
      {
	 gw[l+gl*r][0]=g[l+(gl-glh)+gl*r][0];
	 gw[l+gl*r][1]=-g[l+(gl-glh)+gl*r][1];
      }
      for (l=glh;l<gl;l++)
      {
	 gw[l+gl*r][0]=g[l-glh+gl*r][0];
	 gw[l+gl*r][1]=-g[l-glh+gl*r][1];
      }
   }

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (n=0; n<glh_d_a; n++)
   {
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {

	    fbd=(LTFAT_REAL*)f+2*(L-(glh-n*a)+L*w);
	    for (l=0;l<glh-n*a;l++)
	    {
	       fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }
	    fbd=(LTFAT_REAL*)f-2*(glh-n*a)+2*L*w;
	    for (l=glh-n*a;l<gl;l++)
	    {
	       fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }

	    THE_SUM

	 }	       

      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {

      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=(LTFAT_REAL*)f+2*(n*a-glh+L*w);
	    for (l=0;l<gl;l++)
	    {
	       fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }

	    THE_SUM
	 }

      }

   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=(LTFAT_REAL*)f+2*(n*a-glh+L*w);
	    for (l=0;l<L-n*a+glh;l++)
	    {
	       fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }
	    fbd=(LTFAT_REAL*)f-2*(L-n*a+glh)+2*L*w;
	    for (l=L-n*a+glh;l<gl;l++)
	    {	
	       fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }

	    THE_SUM
	 }
      }
   }

    /* -----------  Clean up ----------------- */   
   ltfat_free(sbuf);    
   ltfat_free(gw);
   ltfat_free(fw);    

   LTFAT_FFTW(destroy_plan)(p_small);
}

/* See the comments on the macro THE_SUM. This macro uses real valued
 * inputs, but produces complex valued output and uses a regular FFT.
 */
#define THE_SUM_R {for (m=0;m<M;m++) \
{ \
   rem = 2*positiverem(m+(n*a-glh), M); \
   sbuf[rem]=0.0; \
   sbuf[rem+1]=0.0; \
   fbd=fw+m; \
   for (k=0;k<gl/M;k++) \
   { \
      sbuf[rem]  += (*fbd); \
      fbd+=M; \
   } \
}   \
\
    LTFAT_FFTW(execute)(p_small);		\
\
coefsum=(LTFAT_REAL*)cout+2*(n*M+r*M*N+w*M*N*R); \
for (m=0;m<M;m++) \
{ \
   coefsum[2*m]   = sbuf[2*m]; \
   coefsum[2*m+1] = sbuf[2*m+1]; \
}} 



void LTFAT_NAME(dgt_fb_r)(const LTFAT_REAL *f, const LTFAT_REAL *g,
	      const int L, const int gl,
	      const int W, const int a, const int M, 
	      LTFAT_COMPLEX *cout)
{
   /*  --------- initial declarations -------------- */

   int k, l, m, n, r, w, rem;

   LTFAT_REAL *gw;

   LTFAT_FFTW(plan) p_small;

   LTFAT_REAL *gb;
   LTFAT_REAL *sbuf,*coefsum, *fw;

   const LTFAT_REAL *fbd;
   

   const int R = 1; 

   /*  ----------- calculation of parameters and plans -------- */
   
   const int N=L/a;

   /* This is a floor operation. */
   const int glh=gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   
   gw   = (LTFAT_REAL*)ltfat_malloc(gl*R*sizeof(LTFAT_REAL));
   fw   = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));
   sbuf = (LTFAT_REAL*)ltfat_malloc(2*M*sizeof(LTFAT_REAL));

   /* Create plan. In-place. */
   p_small = LTFAT_FFTW(plan_dft_1d)(M, (LTFAT_COMPLEX*)sbuf,
				     (LTFAT_COMPLEX*)sbuf,
				     FFTW_FORWARD, FFTW_MEASURE);
   
   /*  ---------- main body ----------- */
  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */
   
   for (r=0;r<R;r++)
   {
      for (l=0;l<glh;l++)
      {
	 gw[l+gl*r]=g[l+(gl-glh)+gl*r];
      }
      for (l=glh;l<gl;l++)
      {
	 gw[l+gl*r]=g[l-glh+gl*r];
      }
   }

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (n=0; n<glh_d_a; n++)
   {
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {

	    fbd=f+L-(glh-n*a)+L*w;
	    for (l=0;l<glh-n*a;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }
	    fbd=f-(glh-n*a)+L*w;
	    for (l=glh-n*a;l<gl;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_R

	 }	       

      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {

      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=f+(n*a-glh+L*w);
	    for (l=0;l<gl;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_R
	 }

      }

   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=f+(n*a-glh+L*w);
	    for (l=0;l<L-n*a+glh;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }
	    fbd=f-(L-n*a+glh)+L*w;
	    for (l=L-n*a+glh;l<gl;l++)
	    {	
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_R
	 }
      }
   }

    /* -----------  Clean up ----------------- */   
   ltfat_free(sbuf);    
   ltfat_free(gw);
   ltfat_free(fw);    

   LTFAT_FFTW(destroy_plan)(p_small);
}

void LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
	      const int L, const int gl,
	      const int W, const int a, const int M, 
	      LTFAT_COMPLEX *cout)
{
   /*  --------- initial declarations -------------- */

   int k, l, m, n, r, w, rem;

   LTFAT_REAL *gw;

   LTFAT_FFTW(plan) p_small;

   LTFAT_REAL *gb;
   LTFAT_REAL *sbuf_in, *coefsum, *fw;
   LTFAT_COMPLEX *sbuf_out;

   const LTFAT_REAL *fbd;
   

   const int R = 1; 

   /*  ----------- calculation of parameters and plans -------- */
   
   const int N=L/a;

   /* These are floor operations. */
   const int glh=gl/2;
   const int M2=M/2+1;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   
   gw   = (LTFAT_REAL*)ltfat_malloc(gl*R*sizeof(LTFAT_REAL));
   fw   = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));
   sbuf_in  = (LTFAT_REAL*)ltfat_malloc(M*sizeof(LTFAT_REAL));
   sbuf_out = (LTFAT_COMPLEX*)ltfat_malloc(M2*sizeof(LTFAT_COMPLEX));

   /* Create plan. In-place. */

   p_small = LTFAT_FFTW(plan_dft_r2c_1d)(M, sbuf_in, sbuf_out, FFTW_MEASURE);
   
   /*  ---------- main body ----------- */
  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */
   
   for (r=0;r<R;r++)
   {
      for (l=0;l<glh;l++)
      {
	 gw[l+gl*r]=g[l+(gl-glh)+gl*r];
      }
      for (l=glh;l<gl;l++)
      {
	 gw[l+gl*r]=g[l-glh+gl*r];
      }
   }

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (n=0; n<glh_d_a; n++)
   {
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {

	    fbd=f+L-(glh-n*a)+L*w;
	    for (l=0;l<glh-n*a;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }
	    fbd=f-(glh-n*a)+L*w;
	    for (l=glh-n*a;l<gl;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_REAL

	 }	       

      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {

      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=f+(n*a-glh+L*w);
	    for (l=0;l<gl;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_REAL
	 }

      }

   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=f+(n*a-glh+L*w);
	    for (l=0;l<L-n*a+glh;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }
	    fbd=f-(L-n*a+glh)+L*w;
	    for (l=L-n*a+glh;l<gl;l++)
	    {	
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_REAL
	 }
      }
   }

    /* -----------  Clean up ----------------- */   
   ltfat_free(sbuf_out);
   ltfat_free(sbuf_in);
   ltfat_free(gw);
   ltfat_free(fw);

   LTFAT_FFTW(destroy_plan)(p_small);
}
