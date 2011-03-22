#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "ltfat.h"
#include "tfutil.c"


void LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
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


   /* % Apply ifft to the coefficients. */
   /* coef=ifft(reshape(coef,M,N*W))*sqrt(M); */

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
	 /* Copy to c-buffer and ifft it */
	 for (int m=0; m<M; m++)
	 {
	    cbuf[m][0]=cin[m+n*M+w*M*N][0];
	    cbuf[m][1]=cin[m+n*M+w*M*N][1];
	 }
	 LTFAT_FFTW(execute)(p_small);

	 const int rem = positiverem(-n*a+glh, M);	 
	 for (int ii=0; ii<gl/M; ii++)
	 {
	    for (int m=0; m<rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],gw[m+ii*M]);
	    }
	    for (int m=0; m<M-rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],gw[m+rem+ii*M]);
	    }
	 }

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
	 /* Copy to c-buffer and ifft it */
	 for (int m=0; m<M; m++)
	 {
	    cbuf[m][0]=cin[m+n*M+w*M*N][0];
	    cbuf[m][1]=cin[m+n*M+w*M*N][1];
	 }
	 LTFAT_FFTW(execute)(p_small);

      	 const int rem = positiverem(-n*a+glh,M);
      	 for (int ii=0; ii<gl/M; ii++)
      	 {
      	    for (int m=0; m<rem; m++)
      	    {
      	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],gw[m+ii*M]);
      	    }
      	    for (int m=0; m<M-rem; m++)
      	    {
      	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],gw[m+rem+ii*M]);
      	    }
      	 }
	 
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

	 /* Copy to c-buffer and ifft it */
	 for (int m=0; m<M; m++)
	 {
	    cbuf[m][0]=cin[m+n*M+w*M*N][0];
	    cbuf[m][1]=cin[m+n*M+w*M*N][1];
	 }
	 LTFAT_FFTW(execute)(p_small);

      	 const int rem = positiverem(-n*a+glh,M);
      	 for (int ii=0; ii<gl/M; ii++)
      	 {
      	    for (int m=0; m<rem; m++)
      	    {
      	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],gw[m+ii*M]);
      	    }
      	    for (int m=0; m<M-rem; m++)
      	    {
      	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],gw[m+rem+ii*M]);
      	    }
      	 }
	 
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
