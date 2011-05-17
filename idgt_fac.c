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
   
   LTFAT_FFTW(plan) p_before, p_after, p_veryend;
   LTFAT_COMPLEX *ff, *cf, *cwork, *cbuf;

   /*  ----------- calculation of parameters and plans -------- */

   const int b=L/M;
   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=b/p;

   h_a=-h_a;

   ff    = (LTFAT_COMPLEX*)ltfat_malloc(d*p*q*W*sizeof(LTFAT_COMPLEX));
   cf    = (LTFAT_COMPLEX*)ltfat_malloc(d*q*q*W*sizeof(LTFAT_COMPLEX));
   cwork = (LTFAT_COMPLEX*)ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEX));
   cbuf  = (LTFAT_COMPLEX*)ltfat_malloc(d*sizeof(LTFAT_COMPLEX));

   /* Scaling constant needed because of FFTWs normalization. */
   const double scalconst = 1.0/((double)d*sqrt((double)M));

   /* Create plans. In-place. */   

   p_after  = LTFAT_FFTW(plan_dft_1d)(d, cbuf, cbuf,
			       FFTW_FORWARD, FFTW_ESTIMATE);

   p_before = LTFAT_FFTW(plan_dft_1d)(d, cbuf, cbuf,
			       FFTW_BACKWARD, FFTW_ESTIMATE);

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

   const int ld4c=M*N;

   /* Leading dimensions of cf */
   const int ld3b=q*q*W;

   /* Leading dimensions of the 4dim array. */
   const int ld2ff=p*q*W;

   for (int r=0;r<c;r++)
   {	


      LTFAT_COMPLEX *cfp=cf;

      for (int w=0;w<W;w++)
      {
	 /* Complete inverse fac of coefficients */
	 for (int l=0;l<q;l++)
	 {
	    for (int u=0;u<q;u++)
	    {	       	       
	       for (int s=0;s<d;s++)	       
	       {	
		  const int rem = r+l*c+positiverem(u+s*q-l*h_a,N)*M+w*ld4c;
		  cbuf[s][0] = cwork[rem][0];
		  cbuf[s][1] = cwork[rem][1];
	       }		    
	       
	       /* Do inverse fft of length d */
	       LTFAT_FFTW(execute)(p_after);

	       for (int s=0;s<d;s++)	       
	       {	
		  cfp[s*ld3b][0] = cbuf[s][0];
		  cfp[s*ld3b][1] = cbuf[s][1];
	       }
	       /* Advance the cf pointer. This is only done in this
		* one place, because the loops are placed such that
		* this pointer will advance linearly through
		* memory. Reordering the loops will break this. */
	       cfp++;
	    }
	 }
      }            

      
      
      /* -------- compute matrix multiplication ---------- */
      
      
      /* Do the matmul  */
      for (int s=0;s<d;s++)
      {	

	 const LTFAT_COMPLEX *gbase = gf+(r+s*c)*p*q;
	 LTFAT_COMPLEX       *fbase = ff+s*p*q*W;
	 const LTFAT_COMPLEX *cbase = (const LTFAT_COMPLEX *)cf+s*q*q*W;
	 
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


      LTFAT_COMPLEX *ffp = ff;
      LTFAT_COMPLEX *fp  = f+r;

      for (int w=0;w<W;w++)
      {
	 for (int l=0;l<q;l++)
	 {
	    for (int k=0;k<p;k++)
	    {
	       for (int s=0;s<d;s++)
	       {		  
		  cbuf[s][0] = ffp[s*ld2ff][0];
		  cbuf[s][1] = ffp[s*ld2ff][1];
	       }
	       	       
	       LTFAT_FFTW(execute)(p_before);
	       
	       for (int s=0;s<d;s++)
	       {		  
		  const int rem = positiverem(k*M+s*p*M-l*h_a*a, L);
		  fp[rem][0] = cbuf[s][0];
		  fp[rem][1] = cbuf[s][1];
	       }
	       
	       /* Advance the ff pointer. This is only done in this
		* one place, because the loops are placed such that
		* this pointer will advance linearly through
		* memory. Reordering the loops will break this. */		  
	       ffp++;
	    }
	 }
	 fp+=L;	 
      }
      fp-=L*W;

      /* ----- Main loop ends here ------------- */
   }           

    /* -----------  Clean up ----------------- */   

   LTFAT_FFTW(destroy_plan)(p_veryend);   
   LTFAT_FFTW(destroy_plan)(p_after);   
   LTFAT_FFTW(destroy_plan)(p_before);   

   ltfat_free(cwork);
   ltfat_free(ff);
   ltfat_free(cf);

   ltfat_free(cbuf);

   
}



LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf, 
			      const int L, const int W,
			      const int a, const int M,
			      LTFAT_REAL *f)
{

   /*  --------- initial declarations -------------- */

   int h_a, h_m;
   
   LTFAT_FFTW(plan) p_before, p_after, p_veryend;
   LTFAT_COMPLEX *ff, *cf, *cbuf;
   LTFAT_REAL *cwork, *sbuf;

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

   ff    = (LTFAT_COMPLEX*)ltfat_malloc(d2*p*q*W*sizeof(LTFAT_COMPLEX));
   cf    = (LTFAT_COMPLEX*)ltfat_malloc(d2*q*q*W*sizeof(LTFAT_COMPLEX));
   cwork = (LTFAT_REAL*)ltfat_malloc(M*N*W*sizeof(LTFAT_REAL));
   cbuf  = (LTFAT_COMPLEX*)ltfat_malloc(d2*sizeof(LTFAT_COMPLEX));
   sbuf  =    (LTFAT_REAL*)ltfat_malloc(d*sizeof(LTFAT_REAL));

   /* Scaling constant needed because of FFTWs normalization. */
   const double scalconst = 1.0/((double)d*sqrt((double)M));

   /* Create plans. In-place. */   
   p_before = LTFAT_FFTW(plan_dft_c2r_1d)(d, cbuf, sbuf, FFTW_ESTIMATE);

   p_after  = LTFAT_FFTW(plan_dft_r2c_1d)(d, sbuf, cbuf, FFTW_ESTIMATE);         
   
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



   const int ld4c=M*N;

   /* Leading dimensions of cf */
   const int ld3b=q*q*W;

   /* Leading dimensions of the 4dim array. */
   const int ld2ff=p*q*W;

   /* -------- Main loop ----------------------------------- */
   for (int r=0;r<c;r++)
   {	

      /* -------- compute coefficient factorization ----------- */

      LTFAT_COMPLEX *cfp=cf;

      for (int w=0;w<W;w++)
      {
	 /* Complete inverse fac of coefficients */
	 for (int l=0;l<q;l++)
	 {
	    for (int u=0;u<q;u++)
	    {	       	       
	       for (int s=0;s<d;s++)	       
	       {	
		  sbuf[s] = cwork[r+l*c+positiverem(u+s*q-l*h_a,N)*M+w*ld4c];
	       }		    
	       
	       /* Do inverse fft of length d */
	       LTFAT_FFTW(execute)(p_after);

	       for (int s=0;s<d2;s++)	       
	       {	
		  cfp[s*ld3b][0] = cbuf[s][0];
		  cfp[s*ld3b][1] = cbuf[s][1];
	       }
	       /* Advance the cf pointer. This is only done in this
		* one place, because the loops are placed such that
		* this pointer will advance linearly through
		* memory. Reordering the loops will break this. */
	       cfp++;
	    }
	 }
      }            


      /* -------- compute matrix multiplication ---------- */
      
      
      /* Do the matmul  */
      for (int s=0;s<d2;s++)
      {	
	 const LTFAT_COMPLEX *gbase = gf+(r+s*c)*p*q;
	 LTFAT_COMPLEX       *fbase = ff+s*p*q*W;
	 const LTFAT_COMPLEX *cbase = (const LTFAT_COMPLEX *)cf+s*q*q*W;
	 
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

      LTFAT_COMPLEX *ffp = ff;
      LTFAT_REAL    *fp  = f+r;

      for (int w=0;w<W;w++)
      {
	 for (int l=0;l<q;l++)
	 {
	    for (int k=0;k<p;k++)
	    {
	       for (int s=0;s<d2;s++)
	       {		  
		  cbuf[s][0] = ffp[s*ld2ff][0];
		  cbuf[s][1] = ffp[s*ld2ff][1];
	       }
	       	       
	       LTFAT_FFTW(execute)(p_before);
	       
	       for (int s=0;s<d;s++)
	       {		  
		  fp[positiverem(k*M+s*p*M-l*h_a*a, L)] = sbuf[s];
	       }
	       
	       /* Advance the ff pointer. This is only done in this
		* one place, because the loops are placed such that
		* this pointer will advance linearly through
		* memory. Reordering the loops will break this. */		  
	       ffp++;
	    }
	 }
	 fp+=L;	 
      }
      fp-=L*W;


      /*  ------- main loop ends -------- */
   }           

    /* -----------  Clean up ----------------- */   

   LTFAT_FFTW(destroy_plan)(p_veryend);   
   LTFAT_FFTW(destroy_plan)(p_after);   
   LTFAT_FFTW(destroy_plan)(p_before);   

   ltfat_free(cwork);
   ltfat_free(ff);
   ltfat_free(cf);

   ltfat_free(cbuf);
   ltfat_free(sbuf);
   
}
