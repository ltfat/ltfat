#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "dgt.h"
#include "tfutil.h"
#include "ltfat_time.h"

/* As opposed to dgt_fb_1, with correct loop ordering and integer optimizations. */

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

   fftw_plan p_small;

   ltfat_complex *fb, *gb;
   double *sbuf,*coefsum, *fbd;

   double st0, st1, st2;
   
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
   fw = (ltfat_complex*)ltfat_malloc(gl*sizeof(ltfat_complex));
   sbuf = (double*)ltfat_malloc(2*M*sizeof(double));

   /* Create plan. In-place. */
   p_small = fftw_plan_dft_1d(M, (ltfat_complex*)sbuf, (ltfat_complex*)sbuf,
			      FFTW_FORWARD, FFTW_MEASURE);
   
   if (dotime)
   {
      st1=ltfat_time();
      printf("DGT_FB_2: Planning phase %f\n",st1-st0);
   }

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
      
      /* Matlab code for this:
       * fpart=[f(L-(glh-n*a)+1:L,:);...
       * f(1:gl-(glh-n*a),:)];
       * 
       * fg=fpart.*gw;
       */
      for (r=0;r<R;r++)
      {
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {

	    fbd=(double*)f+2*(L-(glh-n*a)+L*w);
	    for (l=0;l<glh-n*a;l++)
	    {
	       fw[l][0]=fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[l][1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }
	    fbd=(double*)f-2*(glh-n*a)+2*L*w;
	    for (l=glh-n*a;l<gl;l++)
	    {
	       fw[l][0]=fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[l][1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }

	    /* Do the sum (decimation in frequency, Poisson summation) */      
	    
	    for (m=0;m<M;m++)
	    {
	       rem = 2*positiverem(m+(n*a-glh), M);
	       sbuf[rem]=0.0;
	       sbuf[rem+1]=0.0;
	       fbd=(double*)fw+2*m;
	       for (k=0;k<gl/M;k++)
	       {
		  sbuf[rem]+= fbd[2*k*M];
		  sbuf[rem+1]+= fbd[2*k*M+1];
	       }
	    }
	    
	    fftw_execute(p_small);   	 
	    
	    coefsum=(double*)cout+2*(n*M+r*M*N+w*M*N*R);
	    for (m=0;m<M;m++)
	    {
	       coefsum[2*m] = sbuf[2*m];
	       coefsum[2*m+1] = sbuf[2*m+1];
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
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=(double*)f+2*(n*a-glh+L*w);
	    for (l=0;l<gl;l++)
	    {
	       fw[l][0]=fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[l][1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }

	    /* Do the sum (decimation in frequency, Poisson summation) */      
	    
	    for (m=0;m<M;m++)
	    {
	       rem = 2*positiverem(m+(n*a-glh), M);
	       sbuf[rem]=0.0;
	       sbuf[rem+1]=0.0;

	       fb=fw+m;
	       for (k=0;k<gl/M;k++)
	       {
		  sbuf[rem]+= fb[k*M][0];
		  sbuf[rem+1]+= fb[k*M][1];
	       }
	    }
	    
	    fftw_execute(p_small);   	 
	    
	    coefsum=(double*)cout+2*(n*M+r*M*N+w*M*N*R);
	    for (m=0;m<M;m++)
	    {
	       coefsum[2*m] = sbuf[2*m];
	       coefsum[2*m+1] = sbuf[2*m+1];
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
	 gb=gw+r*gl;
	 for (w=0;w<W;w++)
	 {
	    fbd=(double*)f+2*(n*a-glh+L*w);
	    for (l=0;l<L-n*a+glh;l++)
	    {
	       fw[l][0]=fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[l][1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }
	    fbd=(double*)f-2*(L-n*a+glh)+2*L*w;
	    for (l=L-n*a+glh;l<gl;l++)
	    {	
	       fw[l][0]=fbd[2*l]*gb[l][0]-fbd[2*l+1]*gb[l][1];
	       fw[l][1]=fbd[2*l+1]*gb[l][0]+fbd[2*l]*gb[l][1];
	    }


	    /* Do the sum (decimation in frequency, Poisson summation) */      

	    for (m=0;m<M;m++)
	    {
	       rem = 2*positiverem(m+(n*a-glh), M);
	       sbuf[rem]=0.0;
	       sbuf[rem+1]=0.0;
	       fb=fw+m;
	       for (k=0;k<gl/M;k++)
	       {
		  sbuf[rem]+= fb[k*M][0];
		  sbuf[rem+1]+= fb[k*M][1];
	       }
	    }
	    
	    fftw_execute(p_small);   	 
	 
	    coefsum=(double*)cout+2*(n*M+r*M*N+w*M*N*R);
	    for (m=0;m<M;m++)
	    {
	       coefsum[2*m] = sbuf[2*m];
	       coefsum[2*m+1] = sbuf[2*m+1];
	    }
	 }
      }
   }

   if (dotime)
   {
      st2=ltfat_time();
      printf("DGT_FB_2: Total time %f\n",st2-st0);
   }

    /* -----------  Clean up ----------------- */   
   ltfat_free(sbuf);    
   ltfat_free(gw);
   ltfat_free(fw);    

   fftw_destroy_plan(p_small);
}
