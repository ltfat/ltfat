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

   const int b=L/M;
   const int N=L/a;

   int rem, ep, sp;

   /* This is a floor operation. */
   const int glh=gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   LTFAT_COMPLEX *fw;

   LTFAT_REAL *cbuf = (LTFAT_REAL*)ltfat_malloc(sizeof(LTFAT_COMPLEX));

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
      /* ----- Handle the first boundary using periodic boundary conditions. --- */
      for (int n=0; n<glh_d_a; n++)
      {
	 /* Copy to c-buffer and ifft it */
	 for (int m=0; m<rem; m++)
	 {
	    cbuf[m]=(*cin)[m+n*M+w*M*N];
	 }
	 LTFAT_FFTW(execute)(p_small);

	 rem = positiverem(n*a-glh, M);	 
	 for (int ii=0; ii<gl/M; ii++)
	 {
	    for (int m=0; m<rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],g[m+ii*M]);
	    }
	    for (int m=0; m<M-rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],g[m+rem+ii*M]);
	    }
	 }

	 sp=positiverem(n*a-glh,L);
	 ep=positiverem(n*a-glh+gl-1,L);

	 /* % Add the ff vector to f at position sp. */
	 for (int ii=0;ii<L-sp;ii++}
	 {
	    fw[sp+ii]+=ff[ii];
	 }
	 for (int ii=0; ii<ep+1;ii+)
	 {
	    fw[ii] +=ff[L-sp+ii];
	 }
      }
   

      /* ----- Handle the middle case. --------------------- */
      for (int n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
      {
	 rem=positiverem(-n*a+glh,M);      
	 for (int ii=0; ii<gl/M; ii++)
	 {
	    for (int m=0; m<rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],g[m+ii*M]);
	    }
	    for (int m=0; m<M-rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],g[m+rem+ii*M]);
	    }
	 }
	 
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
	 rem=positiverem(-n*a+glh,M);
	 for (int ii=0; ii<gl/M; ii++)
	 {
	    for (int m=0; m<rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M,cbuf[M-rem+m],g[m+ii*M]);
	    }
	    for (int m=0; m<M-rem; m++)
	    {
	       LTFAT_NAME(complexprod)(ff+m+ii*M+rem,cbuf[m],g[m+rem+ii*M]);
	    }
	 }
	 
	 sp=positiverem(n*a-glh,L);
	 ep=positiverem(n*a-glh+gl-1,L);
	 
	 /* Add the ff vector to f at position sp. */
	 for (int ii=0;ii<L-sp;ii++}
	 {
	    fw[sp+ii]+=ff[ii];
	 }
	 for (int ii=0; ii<ep+1;ii+)
	 {
	    fw[ii] +=ff[L-sp+ii];
	 }
	 
      }
   }
   
   ltfat_free(cbuf);
   ltfat_free(ff);
   ltfat_free(gw);
   
   LTFAT_FFTW(destroy_plan)(p_small);

}
