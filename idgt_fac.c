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
LTFAT_NAME(idgt_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf, 
			  const int L, const int W,
			  const int a, const int M,
			  LTFAT_COMPLEX *f)
{

   /*  --------- initial declarations -------------- */

   int h_a, h_m;
   
   LTFAT_COMPLEX *fbase, *cbase;

   const LTFAT_COMPLEX *gbase;

   div_t domod;

   LTFAT_FFTW(plan) p_before, p_after, p_veryend;
   LTFAT_COMPLEX *ff, *cf, *cwork;

   /*  ----------- calculation of parameters and plans -------- */

   const int b=L/M;
   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   h_a=-h_a;

   ff    = (LTFAT_COMPLEX*)ltfat_malloc(d*p*q*W*sizeof(LTFAT_COMPLEX));
   cf    = (LTFAT_COMPLEX*)ltfat_malloc(q*N*W*sizeof(LTFAT_COMPLEX));
   cwork = (LTFAT_COMPLEX*)ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEX));

   /* Scaling constant needed because of FFTWs normalization. */
   const double scalconst = 1.0/((double)d*sqrt((double)M));

   /* Create plans. In-place. */   
   p_before = LTFAT_FFTW(plan_many_dft)(1, &d, p*q*W,
				 ff, NULL,
				 p*q*W, 1,
				 ff, NULL,
				 p*q*W, 1,
				 FFTW_BACKWARD, FFTW_ESTIMATE);


   p_after = LTFAT_FFTW(plan_many_dft)(1, &d, q*q*W,
				cf, NULL,
				q*q*W, 1,
				cf, NULL,
				q*q*W, 1,
				FFTW_FORWARD, FFTW_ESTIMATE);
   

   /* Create plan. Copy data so we do not overwrite input. Therefore
      it is ok to cast away the constness of cin.*/
   p_veryend = LTFAT_FFTW(plan_many_dft)(1, &M, N*W,
                                  (LTFAT_COMPLEX *)cin, NULL,
                                  1, M,
                                  cwork, NULL,
				  1, M,
				  FFTW_BACKWARD, FFTW_ESTIMATE);

   /* -------- Execute initial IFFT ------------------------ */
   LTFAT_FFTW(execute)(p_veryend);   


   /* -------- Main loop ----------------------------------- */

   for (int r=0;r<c;r++)
   {	


      /* -------- compute coefficient factorization ----------- */

      /* Leading dimensions of the 4dim array. */
      const int ld1=q;
      int ld2=q*q*W;
      
      for (int w=0;w<W;w++)
      {
	 for (int s=0;s<d;s++)
	 {
	    for (int l=0;l<q;l++)
	    {
	       for (int u=0;u<q;u++)
	       {	       
		  /*Add N to make sure it is positive */
		  domod= div(u+s*q-l*h_a+N*M,N);
		  
		  cf[u+(l+q*w)*ld1+s*ld2][0] = cwork[r+l*c+domod.rem*M+w*M*N][0];
		  cf[u+(l+q*w)*ld1+s*ld2][1] = cwork[r+l*c+domod.rem*M+w*M*N][1];
	       }
	    }
	 }
      }           
      
      /* Do fft of length d */
      LTFAT_FFTW(execute)(p_after);
      
      
      /* -------- compute matrix multiplication ---------- */
      
      
      /* Do the matmul  */
      for (int s=0;s<d;s++)
      {	
	 gbase=gf+(r+s*c)*p*q;
	 fbase=ff+s*p*q*W;
	 cbase=cf+s*q*q*W;
	 
	 for (int nm=0;nm<q*W;nm++)
	 {
	    for (int km=0;km<p;km++)
	    {
#ifdef HAVE_COMPLEX_H
	       fbase[km+nm*p]=0.0;
	       for (int mm=0;mm<q;mm++)
	       {
		 fbase[km+nm*p]+=gbase[km+mm*p]*cbase[mm+nm*q];
	       }
	       /* Scale because of FFTWs normalization. */
	       fbase[km+nm*p]=fbase[km+nm*p]*scalconst;
#else
	       fbase[km+nm*p][0]=0.0;
	       fbase[km+nm*p][1]=0.0;
	       for (int mm=0;mm<q;mm++)
	       {
		 fbase[km+nm*p][0]+=gbase[km+mm*p][0]*cbase[mm+nm*q][0]-gbase[km+mm*p][1]*cbase[mm+nm*q][1];
		 fbase[km+nm*p][1]+=gbase[km+mm*p][0]*cbase[mm+nm*q][1]+gbase[km+mm*p][1]*cbase[mm+nm*q][0];
	       }
	       /* Scale because of FFTWs normalization. */
	       fbase[km+nm*p][0]=fbase[km+nm*p][0]*scalconst;
	       fbase[km+nm*p][1]=fbase[km+nm*p][1]*scalconst;
#endif
	    }		  
	 }	      	 
      }

         


      /* ----------- compute inverse signal factorization ---------- */


      /* Do ifft to begin inverse signal factorization.*/
      LTFAT_FFTW(execute)(p_before);

      /* Leading dimensions of the 4dim array. */
      ld2=p*q*W;

      for (int w=0;w<W;w++)
      {
	 for (int s=0;s<d;s++)
	 {
	    for (int l=0;l<q;l++)
	    {
	       for (int k=0;k<p;k++)
	       {
		  /* Add L*M to make sure it is always positive */
		  domod = div(k*M+s*p*M+l*(c-h_m*M)+L*M, L);		  
		  f[r+domod.rem+L*w][0] = ff[k+(l+q*w)*p+s*ld2][0];
		  f[r+domod.rem+L*w][1] = ff[k+(l+q*w)*p+s*ld2][1];
	       }
	    }
	 }
      }


      /* ----- Main loop ends here ------------- */
   }           

    /* -----------  Clean up ----------------- */   

   LTFAT_FFTW(destroy_plan)(p_veryend);   
   LTFAT_FFTW(destroy_plan)(p_after);   
   LTFAT_FFTW(destroy_plan)(p_before);   

   ltfat_free(cwork);
   ltfat_free(ff);
   ltfat_free(cf);
   
}



LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf, 
			      const int L, const int W,
			      const int a, const int M,
			      LTFAT_REAL *f)
{

   /*  --------- initial declarations -------------- */

   int h_a, h_m;
   
   LTFAT_COMPLEX *fbase, *cbase;

   const LTFAT_COMPLEX *gbase;

   div_t domod;

   LTFAT_FFTW(plan) p_before, p_after, p_veryend;
   LTFAT_COMPLEX *ff, *cf;
   LTFAT_REAL *cwork;

   /* This is a floor operation. */
   const int M2= M/2+1;
   
   /*  ----------- calculation of parameters and plans -------- */

   const int b=L/M;
   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   /* This is a floor operation. */
   const int d2= d/2+1;

   h_a=-h_a;

   ff    = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
   cf    = (LTFAT_COMPLEX*)ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEX));
   cwork = (LTFAT_REAL*)ltfat_malloc(M*N*W*sizeof(LTFAT_REAL));

   /* Scaling constant needed because of FFTWs normalization. */
   const double scalconst = 1.0/((double)d*sqrt((double)M));

   /* Create plans. In-place. */   
   p_before = LTFAT_FFTW(plan_many_dft)(1, &d, p*q*W,
				 ff, NULL,
				 p*q*W, 1,
				 ff, NULL,
				 p*q*W, 1,
				 FFTW_BACKWARD, FFTW_ESTIMATE);


   p_after = LTFAT_FFTW(plan_many_dft)(1, &d, q*q*W,
				cf, NULL,
				q*q*W, 1,
				cf, NULL,
				q*q*W, 1,
				FFTW_FORWARD, FFTW_ESTIMATE);
   
   /* Create plan. Copy data so we do not overwrite input. Therefore
      it is ok to cast away the constness of cin. This transform
      destroys its input by default, but the extra flag should prevent
      this. */
   p_veryend = LTFAT_FFTW(plan_many_dft_c2r)(1, &M, N*W,
					    (LTFAT_COMPLEX *)cin, NULL,
					    1, M2,
					    cwork, NULL,
					    1, M,
					    FFTW_ESTIMATE+FFTW_PRESERVE_INPUT);


   /* -------- Execute initial IFFT ------------------------ */
   LTFAT_FFTW(execute)(p_veryend);   


   /* -------- Main loop ----------------------------------- */
   for (int r=0;r<c;r++)
   {	


      /* -------- compute coefficient factorization ----------- */
      
      /* Leading dimensions of the 4dim array. */
      const int ld1=q;
      int ld2=q*q*W;
   
      for (int w=0;w<W;w++)
      {
	 for (int s=0;s<d;s++)
	 {
	    for (int l=0;l<q;l++)
	    {
	       for (int u=0;u<q;u++)
	       {	       
		  /*Add N to make sure it is positive */
		  domod= div(u+s*q-l*h_a+N*M,N);
		  cf[u+(l+q*w)*ld1+s*ld2][0] = cwork[r+l*c+domod.rem*M+w*M*N];
		  cf[u+(l+q*w)*ld1+s*ld2][1] = 0.0;
	       }	       
	    }
	 }
      }           

      
      /* Do fft of length d */
      LTFAT_FFTW(execute)(p_after);


      /* -------- compute matrix multiplication ---------- */
      
      
      /* Do the matmul  */
      for (int s=0;s<d;s++)
      {	
	 gbase=gf+(r+s*c)*p*q;
	 fbase=ff+s*p*q*W;
	 cbase=cf+s*q*q*W;
	 
	 for (int nm=0;nm<q*W;nm++)
	 {
	    for (int km=0;km<p;km++)
	    {
#ifdef HAVE_COMPLEX_H
	       fbase[km+nm*p]=0.0;
	       for (int mm=0;mm<q;mm++)
	       {
		 fbase[km+nm*p]+=gbase[km+mm*p]*cbase[mm+nm*q];
	       }
	       /* Scale because of FFTWs normalization. */
	       fbase[km+nm*p]=fbase[km+nm*p]*scalconst;
#else
	       fbase[km+nm*p][0]=0.0;
	       fbase[km+nm*p][1]=0.0;
	       for (int mm=0;mm<q;mm++)
	       {
		 fbase[km+nm*p][0]+=gbase[km+mm*p][0]*cbase[mm+nm*q][0]-gbase[km+mm*p][1]*cbase[mm+nm*q][1];
		 fbase[km+nm*p][1]+=gbase[km+mm*p][0]*cbase[mm+nm*q][1]+gbase[km+mm*p][1]*cbase[mm+nm*q][0];
	       }
	       /* Scale because of FFTWs normalization. */
	       fbase[km+nm*p][0]=fbase[km+nm*p][0]*scalconst;
	       fbase[km+nm*p][1]=fbase[km+nm*p][1]*scalconst;
#endif
	    }		  
	 }	      	 
      }        


      /* ----------- compute inverse signal factorization ---------- */


      /* Do ifft to begin inverse signal factorization.*/
      LTFAT_FFTW(execute)(p_before);

      /* Leading dimensions of the 4dim array. */
      ld2=p*q*W;

      for (int w=0;w<W;w++)
      {
	 for (int s=0;s<d;s++)
	 {
	    for (int l=0;l<q;l++)
	    {
	       for (int k=0;k<p;k++)
	       {
		  /* Add L*M to make sure it is always positive */
		  domod = div(k*M+s*p*M+l*(c-h_m*M)+L*M, L);		  
		  f[r+domod.rem+L*w] = ff[k+(l+q*w)*p+s*ld2][0];
	       }
	    }
	 }
      }

      /*  ------- main loop ends -------- */
   }           

    /* -----------  Clean up ----------------- */   

   LTFAT_FFTW(destroy_plan)(p_veryend);   
   LTFAT_FFTW(destroy_plan)(p_after);   
   LTFAT_FFTW(destroy_plan)(p_before);   

   ltfat_free(cwork);
   ltfat_free(ff);
   ltfat_free(cf);
   
}
