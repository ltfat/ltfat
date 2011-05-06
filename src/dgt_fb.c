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
for (int m=0;m<M;m++) \
{ \
   rem = 2*positiverem(m+(n*a-glh), M); \
   sbuf[rem]=0.0; \
   sbuf[rem+1]=0.0; \
   fbd=fw+2*m; \
   for (int k=0;k<gl/M;k++) \
   { \
      sbuf[rem]  += fbd[0]; \
      sbuf[rem+1]+= fbd[1]; \
      fbd+=2*M; \
   } \
} \
\
 LTFAT_FFTW(execute)(plan.p_small);			\
\
coefsum=(LTFAT_REAL*)cout+2*(n*M+w*M*N); \
for (int m=0;m<M;m++) \
{ \
   coefsum[2*m] = sbuf[2*m]; \
   coefsum[2*m+1] = sbuf[2*m+1]; \
}} 

#define THE_SUM_REAL { \
for (int m=0;m<M;m++) \
{ \
   rem = positiverem(m+(n*a-glh), M); \
   sbuf[rem]=0.0; \
   fbd=fw+m; \
   for (int k=0;k<gl/M;k++) \
   { \
     sbuf[rem]+=(*fbd);			\
      fbd+=M; \
   } \
} \
\
 LTFAT_FFTW(execute)(plan.p_small);			\
\
coefsum=(LTFAT_REAL*)cout+2*(n*M2+w*M2*N); \
for (int m=0;m<M2;m++) \
{ \
   coefsum[2*m]   = cbuf[m][0]; \
   coefsum[2*m+1] = cbuf[m][1]; \
}}


LTFAT_EXTERN LTFAT_NAME(dgt_fb_plan)
LTFAT_NAME(dgt_fb_init)(const LTFAT_COMPLEX *g,
   const int gl, const int a, const int M, 
   unsigned flags)
{
   LTFAT_NAME(dgt_fb_plan) plan;
  
   plan.a=a;
   plan.M=M;
   plan.gl=gl;

   plan.gw  = (LTFAT_COMPLEX*)ltfat_malloc(gl*sizeof(LTFAT_COMPLEX));

   plan.fw  = (LTFAT_REAL*)ltfat_malloc(2*gl*sizeof(LTFAT_REAL));

   plan.sbuf = (LTFAT_REAL*)ltfat_malloc(M*sizeof(LTFAT_COMPLEX));

   plan.p_small = LTFAT_FFTW(plan_dft_1d)(M, (LTFAT_COMPLEX*)plan.sbuf, (LTFAT_COMPLEX*)plan.sbuf,
			      FFTW_FORWARD, flags);

   /* This is a floor operation. */
   const int glh=gl/2;

  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */
   
   for (int l=0;l<glh;l++)
   {
      plan.gw[l][0]=g[l+(gl-glh)][0];
      plan.gw[l][1]=-g[l+(gl-glh)][1];
   }
   for (int l=glh;l<gl;l++)
   {
      plan.gw[l][0]=g[l-glh][0];
      plan.gw[l][1]=-g[l-glh][1];
   }

   return (plan);
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_done)(LTFAT_NAME(dgt_fb_plan) plan)
{

   ltfat_free(plan.sbuf);    
   ltfat_free(plan.gw);
   ltfat_free(plan.fw);    

   LTFAT_FFTW(destroy_plan)(plan.p_small);

}



LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_execute)(LTFAT_NAME(dgt_fb_plan) plan, const LTFAT_COMPLEX *f, 
	      const int L, const int W,  LTFAT_COMPLEX *cout)
{
   /*  --------- initial declarations -------------- */

   const int a=plan.a;
   const int M=plan.M;
   const int N=L/a;
 
   const int gl=plan.gl;
   LTFAT_REAL *sbuf=plan.sbuf;
   LTFAT_REAL *fw=plan.fw;

   /* This is a floor operation. */
   const int glh=plan.gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));

   int rem;

   LTFAT_COMPLEX *gb;
   LTFAT_REAL *coefsum, *fbd;
   
      
   /*  ---------- main body ----------- */

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (int n=0; n<glh_d_a; n++)
   {
      gb=plan.gw;
      for (int w=0;w<W;w++)
      {
	 
	 fbd=(LTFAT_REAL*)f+2*(L-(glh-n*a)+L*w);
	 for (int l=0;l<glh-n*a;l++)
	 {
	    fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	    fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	 }
	 fbd=(LTFAT_REAL*)f-2*(glh-n*a)+2*L*w;
	 for (int l=glh-n*a;l<gl;l++)
	 {
	    fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	    fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	 }
	 
	 THE_SUM

      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {
      gb=plan.gw;
      for (int w=0;w<W;w++)
      {
	 fbd=(LTFAT_REAL*)f+2*(n*a-glh+L*w);
	 for (int l=0;l<gl;l++)
	 {
	    fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	    fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	 }
	 
	 THE_SUM
      }

   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (int n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      gb=plan.gw;
      for (int w=0;w<W;w++)
      {
	 fbd=(LTFAT_REAL*)f+2*(n*a-glh+L*w);
	 for (int l=0;l<L-n*a+glh;l++)
	 {
	    fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	    fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	 }
	 fbd=(LTFAT_REAL*)f-2*(L-n*a+glh)+2*L*w;
	 for (int l=L-n*a+glh;l<gl;l++)
	 {	
	    fw[2*l]  =fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	    fw[2*l+1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	 }
	 
	 THE_SUM
      }
   }

}

/* See the comments on the macro THE_SUM. This macro uses real valued
 * inputs, but produces complex valued output and uses a regular FFT.
 */
#define THE_SUM_R {for (int m=0;m<M;m++) \
{ \
   rem = 2*positiverem(m+(n*a-glh), M); \
   sbuf[rem]=0.0; \
   sbuf[rem+1]=0.0; \
   fbd=fw+m; \
   for (int k=0;k<gl/M;k++) \
   { \
      sbuf[rem]  += (*fbd); \
      fbd+=M; \
   } \
}   \
\
    LTFAT_FFTW(execute)(p_small);		\
\
coefsum=(LTFAT_REAL*)cout+2*(n*M+r*M*N+w*M*N*R); \
for (int m=0;m<M;m++) \
{ \
   coefsum[2*m]   = sbuf[2*m]; \
   coefsum[2*m+1] = sbuf[2*m+1]; \
}} 


LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_r)(const LTFAT_REAL *f, const LTFAT_REAL *g,
	      const int L, const int gl,
	      const int W, const int a, const int M, 
	      LTFAT_COMPLEX *cout)
{
   /*  --------- initial declarations -------------- */

   int r, rem;

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
      for (int l=0;l<glh;l++)
      {
	 gw[l+gl*r]=g[l+(gl-glh)+gl*r];
      }
      for (int l=glh;l<gl;l++)
      {
	 gw[l+gl*r]=g[l-glh+gl*r];
      }
   }

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (int n=0; n<glh_d_a; n++)
   {
      for (int r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (int w=0;w<W;w++)
	 {

	    fbd=f+L-(glh-n*a)+L*w;
	    for (int l=0;l<glh-n*a;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }
	    fbd=f-(glh-n*a)+L*w;
	    for (int l=glh-n*a;l<gl;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_R

	 }	       

      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {

      for (int r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (int w=0;w<W;w++)
	 {
	    fbd=f+(n*a-glh+L*w);
	    for (int l=0;l<gl;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }

	    THE_SUM_R
	 }

      }

   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (int n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      for (int r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (int w=0;w<W;w++)
	 {
	    fbd=f+(n*a-glh+L*w);
	    for (int l=0;l<L-n*a+glh;l++)
	    {
	       fw[l]  =fbd[l]*gb[l];
	    }
	    fbd=f-(L-n*a+glh)+L*w;
	    for (int l=L-n*a+glh;l<gl;l++)
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








LTFAT_EXTERN LTFAT_NAME(dgtreal_fb_plan)
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL *g,
   const int gl, const int a, const int M, 
   unsigned flags)
{
   LTFAT_NAME(dgtreal_fb_plan) plan;
  
   plan.a=a;
   plan.M=M;
   const int M2=M/2+1;
   plan.gl=gl;



   plan.gw   = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));

   plan.fw   = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));

   plan.sbuf = (LTFAT_REAL*)ltfat_malloc(M*sizeof(LTFAT_REAL));

   plan.cbuf = (LTFAT_COMPLEX*)ltfat_malloc(M2*sizeof(LTFAT_COMPLEX));

   plan.p_small = LTFAT_FFTW(plan_dft_r2c_1d)(M, plan.sbuf, plan.cbuf, flags);

   /* This is a floor operation. */
   const int glh=gl/2;

  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */

   for (int l=0;l<glh;l++)
   {
      plan.gw[l]=g[l+(gl-glh)];
   }
   for (int l=glh;l<gl;l++)
   {
      plan.gw[l]=g[l-glh];
   }   

   return (plan);
}

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan) plan)
{

   ltfat_free(plan.sbuf);
   ltfat_free(plan.cbuf);
   ltfat_free(plan.gw);
   ltfat_free(plan.fw);    

   LTFAT_FFTW(destroy_plan)(plan.p_small);

}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_execute)(const LTFAT_NAME(dgtreal_fb_plan) plan, const LTFAT_REAL *f,
			       const int L, const int W, 
			       LTFAT_COMPLEX *cout)
{
   /*  --------- initial declarations -------------- */

   const int a=plan.a;
   const int M=plan.M;
   const int N=L/a;
 
   const int gl=plan.gl;
   LTFAT_REAL    *sbuf=plan.sbuf;
   LTFAT_COMPLEX *cbuf=plan.cbuf;
   LTFAT_REAL *fw=plan.fw;

   /* These are floor operations. */
   const int glh=plan.gl/2;
   const int M2=M/2+1;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));

   int rem;

   LTFAT_REAL *gb;
   LTFAT_REAL *coefsum;

   const LTFAT_REAL *fbd;

   /*  ---------- main body ----------- */

   /*----- Handle the first boundary using periodic boundary conditions.*/
   for (int n=0; n<glh_d_a; n++)
   {
      gb=plan.gw;
      for (int w=0;w<W;w++)
      {
	 
	 fbd=f+L-(glh-n*a)+L*w;
	 for (int l=0;l<glh-n*a;l++)
	 {
	    fw[l]  =fbd[l]*gb[l];
	 }
	 fbd=f-(glh-n*a)+L*w;
	 for (int l=glh-n*a;l<gl;l++)
	 {
	    fw[l]  =fbd[l]*gb[l];
	 }
	 
	 THE_SUM_REAL

      }	     
   }
   
   /* ----- Handle the middle case. --------------------- */
   for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
   {

      gb=plan.gw;
      for (int w=0;w<W;w++)
      {
	 fbd=f+(n*a-glh+L*w);
	 for (int l=0;l<gl;l++)
	 {
	    fw[l]  =fbd[l]*gb[l];
	 }
	 
	 THE_SUM_REAL
      }

   }
   
   /* Handle the last boundary using periodic boundary conditions. */   
   for (int n=(L-(gl+1)/2)/a+1; n<N; n++)
   {
      gb=plan.gw;
      for (int w=0;w<W;w++)
      {
	 fbd=f+(n*a-glh+L*w);
	 for (int l=0;l<L-n*a+glh;l++)
	 {
	    fw[l]  =fbd[l]*gb[l];
	 }
	 fbd=f-(L-n*a+glh)+L*w;
	 for (int l=L-n*a+glh;l<gl;l++)
	 {	
	    fw[l]  =fbd[l]*gb[l];
	 }
	 
         THE_SUM_REAL
      }
   }

}
