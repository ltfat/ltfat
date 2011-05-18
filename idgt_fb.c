#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"

#define THE_SUM { \
	 /* Copy to c-buffer and ifft it */ \
	 for (int m=0; m<M; m++) \
	 { \
	    cbuf[m][0]=cin[m+n*M+w*M*N][0]; \
	    cbuf[m][1]=cin[m+n*M+w*M*N][1]; \
	 } \
	 LTFAT_FFTW(execute)(p_small); \
\
	 const int rem = positiverem(-n*a+glh, M); \
	 for (int ii=0; ii<gl/M; ii++) \
	 { \
	    for (int m=0; m<rem; m++) \
	    { \
	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],gw[m+ii*M]); \
	    } \
	    for (int m=0; m<M-rem; m++) \
	    { \
	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],gw[m+rem+ii*M]); \
	    } \
	 } \
}

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
			 const int L, const int gl, const int W,
			 const int a, const int M,
			 LTFAT_COMPLEX *f)

{ 
  /*  --------- initial declarations -------------- */

   const int N=L/a;

   int ep, sp;

   /* This is a floor operation. */
   const int glh=gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   LTFAT_COMPLEX *fw;

   LTFAT_COMPLEX *cbuf = (LTFAT_COMPLEX*)ltfat_malloc(M*sizeof(LTFAT_COMPLEX));

   /* Create plan. In-place. */
   LTFAT_FFTW(plan) p_small = LTFAT_FFTW(plan_dft_1d)(M, cbuf, cbuf, FFTW_BACKWARD, FFTW_MEASURE);

   /* % The fftshift actually makes some things easier. */
   LTFAT_COMPLEX *gw  = (LTFAT_COMPLEX*)ltfat_malloc(gl*sizeof(LTFAT_COMPLEX));
   for (int l=0;l<glh;l++)
   {
      gw[l][0] = g[l+(gl-glh)][0];
      gw[l][1] = g[l+(gl-glh)][1];
   }
   for (int l=glh;l<gl;l++)
   {
      gw[l][0] = g[l-glh][0];
      gw[l][1] = g[l-glh][1];
   }
   
   LTFAT_COMPLEX *ff  = (LTFAT_COMPLEX*)ltfat_malloc(gl*sizeof(LTFAT_COMPLEX));
   
   for (int w=0; w<W; w++)
   {
      fw=f+w*L;
      for (int l=0;l<L;l++)
      {
	 fw[l][0]=0.0;
	 fw[l][1]=0.0;
      }
      /* ----- Handle the first boundary using periodic boundary conditions. --- */
      for (int n=0; n<glh_d_a; n++)
      {

	 THE_SUM;

	 sp=positiverem(n*a-glh,L);
	 ep=positiverem(n*a-glh+gl-1,L);

	 /* % Add the ff vector to f at position sp. */
	 for (int ii=0;ii<L-sp;ii++)
	 {
	    fw[sp+ii][0]+=ff[ii][0];
	    fw[sp+ii][1]+=ff[ii][1];
	 }
	 for (int ii=0; ii<ep+1;ii++)
	 {
	    fw[ii][0] +=ff[L-sp+ii][0];
	    fw[ii][1] +=ff[L-sp+ii][1];
	 }
      }
   

      /* ----- Handle the middle case. --------------------- */
      for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
      {

	 THE_SUM;

      	 sp=positiverem(n*a-glh,L);
      	 ep=positiverem(n*a-glh+gl-1,L);
	 
      	 /* Add the ff vector to f at position sp. */
      	 for (int ii=0;ii<ep-sp+1;ii++)
      	 {
      	    fw[ii+sp][0] += ff[ii][0];
      	    fw[ii+sp][1] += ff[ii][1];
      	 }
      }
      
      /* Handle the last boundary using periodic boundary conditions. */
      for (int n=(L-(gl+1)/2)/a+1; n<N; n++)
      {

	 THE_SUM;
	 
      	 sp=positiverem(n*a-glh,L);
      	 ep=positiverem(n*a-glh+gl-1,L);
	 
      	 /* Add the ff vector to f at position sp. */
      	 for (int ii=0;ii<L-sp;ii++)
      	 {
      	    fw[sp+ii][0]+=ff[ii][0];
      	    fw[sp+ii][1]+=ff[ii][1];
      	 }
      	 for (int ii=0; ii<ep+1;ii++)
      	 {
      	    fw[ii][0] +=ff[L-sp+ii][0];
      	    fw[ii][1] +=ff[L-sp+ii][1];
      	 }
	 
      }
   }
   
   ltfat_free(cbuf);
   ltfat_free(ff);
   ltfat_free(gw);
   
   LTFAT_FFTW(destroy_plan)(p_small);

}

/* ------------------- IDGT_FB_R --------------------- */

#define THE_SUM_R { \
	 /* Copy to c-buffer and ifft it */ \
	 for (int m=0; m<M; m++) \
	 { \
	    cbuf[m][0]=cin[m+n*M+w*M*N][0]; \
	    cbuf[m][1]=cin[m+n*M+w*M*N][1]; \
	 } \
	 LTFAT_FFTW(execute)(p_small); \
\
	 const int rem = positiverem(-n*a+glh, M); \
	 for (int ii=0; ii<gl/M; ii++) \
	 { \
	    for (int m=0; m<rem; m++) \
	    { \
	       ff[m+ii*M][0] = cbuf[M-rem+m][0]*gw[m+ii*M]; \
	       ff[m+ii*M][1] = cbuf[M-rem+m][1]*gw[m+ii*M]; \
	    } \
	    for (int m=0; m<M-rem; m++) \
	    { \
	       ff[m+ii*M+rem][0] = cbuf[m][0]*gw[m+rem+ii*M]; \
	       ff[m+ii*M+rem][1] = cbuf[m][1]*gw[m+rem+ii*M]; \
	    } \
	 } \
}

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb_r)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
			 const int L, const int gl, const int W,
			 const int a, const int M,
			 LTFAT_COMPLEX *f)

{ 
  /*  --------- initial declarations -------------- */

   const int N=L/a;

   int ep, sp;

   /* This is a floor operation. */
   const int glh=gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   LTFAT_COMPLEX *fw;

   LTFAT_COMPLEX *cbuf = (LTFAT_COMPLEX*)ltfat_malloc(M*sizeof(LTFAT_COMPLEX));

   /* Create plan. In-place. */
   LTFAT_FFTW(plan) p_small = LTFAT_FFTW(plan_dft_1d)(M, cbuf, cbuf, FFTW_BACKWARD, FFTW_MEASURE);

   /* % The fftshift actually makes some things easier. */
   LTFAT_REAL *gw  = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));
   for (int l=0;l<glh;l++)
   {
      gw[l] = g[l+(gl-glh)];
   }
   for (int l=glh;l<gl;l++)
   {
      gw[l] = g[l-glh];
   }
   
   LTFAT_COMPLEX *ff  = (LTFAT_COMPLEX*)ltfat_malloc(gl*sizeof(LTFAT_COMPLEX));
   
   for (int w=0; w<W; w++)
   {
      fw=f+w*L;
      for (int l=0;l<L;l++)
      {
	 fw[l][0]=0.0;
	 fw[l][1]=0.0;
      }
      /* ----- Handle the first boundary using periodic boundary conditions. --- */
      for (int n=0; n<glh_d_a; n++)
      {

	 THE_SUM_R;

	 sp=positiverem(n*a-glh,L);
	 ep=positiverem(n*a-glh+gl-1,L);

	 /* % Add the ff vector to f at position sp. */
	 for (int ii=0;ii<L-sp;ii++)
	 {
	    fw[sp+ii][0]+=ff[ii][0];
	    fw[sp+ii][1]+=ff[ii][1];
	 }
	 for (int ii=0; ii<ep+1;ii++)
	 {
	    fw[ii][0] +=ff[L-sp+ii][0];
	    fw[ii][1] +=ff[L-sp+ii][1];
	 }
      }
   

      /* ----- Handle the middle case. --------------------- */
      for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
      {

	 THE_SUM_R;

      	 sp=positiverem(n*a-glh,L);
      	 ep=positiverem(n*a-glh+gl-1,L);
	 
      	 /* Add the ff vector to f at position sp. */
      	 for (int ii=0;ii<ep-sp+1;ii++)
      	 {
      	    fw[ii+sp][0] += ff[ii][0];
      	    fw[ii+sp][1] += ff[ii][1];
      	 }
      }
      
      /* Handle the last boundary using periodic boundary conditions. */
      for (int n=(L-(gl+1)/2)/a+1; n<N; n++)
      {

	 THE_SUM_R;
	 
      	 sp=positiverem(n*a-glh,L);
      	 ep=positiverem(n*a-glh+gl-1,L);
	 
      	 /* Add the ff vector to f at position sp. */
      	 for (int ii=0;ii<L-sp;ii++)
      	 {
      	    fw[sp+ii][0]+=ff[ii][0];
      	    fw[sp+ii][1]+=ff[ii][1];
      	 }
      	 for (int ii=0; ii<ep+1;ii++)
      	 {
      	    fw[ii][0] +=ff[L-sp+ii][0];
      	    fw[ii][1] +=ff[L-sp+ii][1];
      	 }
	 
      }
   }
   
   ltfat_free(cbuf);
   ltfat_free(ff);
   ltfat_free(gw);
   
   LTFAT_FFTW(destroy_plan)(p_small);

}




/* ------------------- IDGTREAL ---------------------- */

#define THE_SUM_REAL { \
	 /* Copy to c-buffer and ifft it */ \
	 for (int m=0; m<M2; m++) \
	 { \
	    cbuf[m][0]=cin[m+n*M2+w*M2*N][0]; \
	    cbuf[m][1]=cin[m+n*M2+w*M2*N][1]; \
	 } \
	 LTFAT_FFTW(execute)(p_small); \
\
	 const int rem = positiverem(-n*a+glh, M); \
	 for (int ii=0; ii<gl/M; ii++) \
	 { \
	    for (int m=0; m<rem; m++) \
	    { \
	       ff[m+ii*M] = crbuf[M-rem+m]*gw[m+ii*M]; \
	    } \
	    for (int m=0; m<M-rem; m++) \
	    { \
	       ff[m+ii*M+rem] = crbuf[m]*gw[m+rem+ii*M]; \
	    } \
	 } \
}

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
			 const int L, const int gl, const int W,
			 const int a, const int M,
			 LTFAT_REAL *f)

{ 
  /*  --------- initial declarations -------------- */

   const int N=L/a;

   /* This is a floor operation. */
   const int M2= M/2+1;

   int ep, sp;

   /* This is a floor operation. */
   const int glh=gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   LTFAT_REAL *fw;

   LTFAT_COMPLEX *cbuf  = (LTFAT_COMPLEX*)ltfat_malloc(M2*sizeof(LTFAT_COMPLEX));
   LTFAT_REAL    *crbuf =    (LTFAT_REAL*)ltfat_malloc( M*sizeof(LTFAT_REAL));

   /* Create plan. In-place. */
   LTFAT_FFTW(plan) p_small = LTFAT_FFTW(plan_dft_c2r_1d)(M, cbuf, crbuf, FFTW_MEASURE);

   /* % The fftshift actually makes some things easier. */
   LTFAT_REAL *gw  = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));
   for (int l=0;l<glh;l++)
   {
      gw[l] = g[l+(gl-glh)];
   }
   for (int l=glh;l<gl;l++)
   {
      gw[l] = g[l-glh];
   }
   
   LTFAT_REAL *ff  = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));
   
   for (int w=0; w<W; w++)
   {
      fw=f+w*L;
      for (int l=0;l<L;l++)
      {
	 fw[l]=0.0;
      }
      /* ----- Handle the first boundary using periodic boundary conditions. --- */
      for (int n=0; n<glh_d_a; n++)
      {

	 THE_SUM_REAL;

	 sp=positiverem(n*a-glh,L);
	 ep=positiverem(n*a-glh+gl-1,L);

	 /* % Add the ff vector to f at position sp. */
	 for (int ii=0;ii<L-sp;ii++)
	 {
	    fw[sp+ii]+=ff[ii];
	 }
	 for (int ii=0; ii<ep+1;ii++)
	 {
	    fw[ii] +=ff[L-sp+ii];
	 }
      }
   

      /* ----- Handle the middle case. --------------------- */
      for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
      {

	 THE_SUM_REAL;

      	 sp=positiverem(n*a-glh,L);
      	 ep=positiverem(n*a-glh+gl-1,L);
	 
      	 /* Add the ff vector to f at position sp. */
      	 for (int ii=0;ii<ep-sp+1;ii++)
      	 {
      	    fw[ii+sp] += ff[ii];
      	 }
      }
      
      /* Handle the last boundary using periodic boundary conditions. */
      for (int n=(L-(gl+1)/2)/a+1; n<N; n++)
      {

	 THE_SUM_REAL;
	 
      	 sp=positiverem(n*a-glh,L);
      	 ep=positiverem(n*a-glh+gl-1,L);
	 
      	 /* Add the ff vector to f at position sp. */
      	 for (int ii=0;ii<L-sp;ii++)
      	 {
      	    fw[sp+ii]+=ff[ii];
      	 }
      	 for (int ii=0; ii<ep+1;ii++)
      	 {
      	    fw[ii] +=ff[L-sp+ii];
      	 }
	 
      }
   }
   
   ltfat_free(cbuf);
   ltfat_free(crbuf);
   ltfat_free(ff);
   ltfat_free(gw);
   
   LTFAT_FFTW(destroy_plan)(p_small);

}
