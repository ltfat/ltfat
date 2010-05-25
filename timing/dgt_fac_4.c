#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "dgt.h"
#include "tfutil.h"
#include "ltfat_time.h"

/*  This routine computes the DGT factorization using strided FFTs so
    the memory layout is optimized for the matrix product. Compared to
    dgt_fac_1, it moves the r-loop to be the outermost loop to
    conserve memory and hopefully use the cache hierachy better */

void dgt_fac(ltfat_complex *f, ltfat_complex *gf, const int L, const int W,
	     const int R, const int a, const int M, ltfat_complex *cout, int dotime)
{

   /*  --------- initial declarations -------------- */

   int b, N, c, d, p, q, h_a, h_m;
   
   ltfat_complex *gbase, *fbase, *cbase;

   int l, k, r, s, u, w, rw, nm, mm, km;
   int ld2a, ld1b, ld3b;
   int rem;

   fftw_plan p_before, p_after, p_veryend;
   ltfat_complex *ff, *cf;

   double scalconst;
   
   double st0, st1, st6, st7;

   /*  ----------- calculation of parameters and plans -------- */
   
   if (dotime)
   {
      st0=ltfat_time();
   }

   b=L/M;
   N=L/a;
   
   c=gcd(a, M,&h_a, &h_m);
   p=a/c;
   q=M/c;
   d=b/p;

   h_a=-h_a;

   /* Scaling constant needed because of FFTWs normalization. */
   scalconst=1.0/((double)d*sqrt((double)M));

   /*printf("%i %i %i %i %i\n",c,d,p,q,W);*/

   ff = (ltfat_complex*)ltfat_malloc(d*p*q*W*sizeof(ltfat_complex));
   cf = (ltfat_complex*)ltfat_malloc(d*q*q*W*R*sizeof(ltfat_complex));

   /* Create plans. In-place. */
   
   p_before = fftw_plan_many_dft(1, &d, p*q*W,
				 ff, NULL,
				 p*q*W, 1,
				 ff, NULL,
				 p*q*W, 1,
				 FFTW_FORWARD, FFTW_OPTITYPE);


   p_after = fftw_plan_many_dft(1, &d, q*q*W*R,
				cf, NULL,
				q*q*W*R, 1,
				cf, NULL,
				q*q*W*R, 1,
				FFTW_BACKWARD, FFTW_OPTITYPE);
   
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
      printf("DGT_FAC_4: Planning phase %f\n",st1-st0);
   }


   /* Leading dimensions of the 4dim array. */
   ld2a=p*q*W;

   /* Leading dimensions of cf */
   ld1b=q*R;
   ld3b=q*R*q*W;
   
   /* --------- main loop begins here ------------------- */
   for (r=0;r<c;r++)
   {


   /*  ---------- compute signal factorization ----------- */

      for (s=0;s<d;s++)
      {		  
	 for (w=0;w<W;w++)
	 {
	    for (l=0;l<q;l++)
	    {
	       for (k=0;k<p;k++)
	       {
		  rem = positiverem(k*M+s*p*M-l*h_a*a, L);
	       
#ifdef HAVE_COMPLEX_H
		  ff[k+(l+q*w)*p+s*ld2a]=f[r+rem+L*w]*scalconst;
#else
		  ff[k+(l+q*w)*p+s*ld2a][0]=f[r+rem+L*w][0]*scalconst;
		  ff[k+(l+q*w)*p+s*ld2a][1]=f[r+rem+L*w][1]*scalconst;
#endif
	       }
	    }
	 }
      }

   
      /* Do fft to complete signal factorization.*/
      fftw_execute(p_before);


      /* ----------- compute matrix multiplication ----------- */

      
      /* Do the matmul  */
      for (s=0;s<d;s++)
      {
	 gbase=gf+(r+s*c)*p*q*R;
	 fbase=ff+s*p*q*W;
	 cbase=cf+s*q*q*W*R;
	 
	 for (nm=0;nm<q*W;nm++)
	 {
	    for (mm=0;mm<q*R;mm++)
	    {
#ifdef HAVE_COMPLEX_H
	       cbase[mm+nm*q*R]=0.0;
	       for (km=0;km<p;km++)
	       {
		 cbase[mm+nm*q*R]+=conj(gbase[km+mm*p])*fbase[km+nm*p];
	       }
#else
	       cbase[mm+nm*q*R][0]=0.0;
	       cbase[mm+nm*q*R][1]=0.0;
	       for (km=0;km<p;km++)
	       {
		  cbase[mm+nm*q*R][0]+=gbase[km+mm*p][0]*fbase[km+nm*p][0]+gbase[km+mm*p][1]*fbase[km+nm*p][1];
		  cbase[mm+nm*q*R][1]+=gbase[km+mm*p][0]*fbase[km+nm*p][1]-gbase[km+mm*p][1]*fbase[km+nm*p][0];
	       }
#endif
	    }		  
	 }	      	 
      }


      /*  -------  compute inverse coefficient factorization ------- */
      
      /* Do inverse fft of length d */
      fftw_execute(p_after);

               
      for (w=0;w<W;w++)
      {
	 /* Complete inverse fac of coefficients */
	 for (rw=0;rw<R;rw++)
	 {
	    for (s=0;s<d;s++)	       
	    {	
	       for (u=0;u<q;u++)
	       {	       
		  
		 for (l=0;l<q;l++)
		 {		    
		     rem= positiverem(u+s*q-l*h_a,N)*M;
#ifdef HAVE_COMPLEX_H	  
		     cout[r+l*c+rem+rw*M*N+w*M*N*R]=cf[u+rw*q+(l+q*w)*ld1b+s*ld3b];
#else
		     cout[r+l*c+rem+rw*M*N+w*M*N*R][0]=cf[u+rw*q+(l+q*w)*ld1b+s*ld3b][0];
		     cout[r+l*c+rem+rw*M*N+w*M*N*R][1]=cf[u+rw*q+(l+q*w)*ld1b+s*ld3b][1];
#endif
		  }
	       }
	    }
	 }
      }      

      
      /* ----------- Main loop ends here ------------------------ */
   }     

   if (dotime)
   {
      st6=ltfat_time();
      printf("DGT_FAC_4: Main loop done %f\n",st6-st1);
   }

   /* FFT to modulate the coefficients. */
   fftw_execute(p_veryend);   

   if (dotime)
   {
      st7=ltfat_time();
      printf("DGT_FAC_4: Final FFT %f\n",st7-st6);
      printf("DGT_FAC_4: Total time %f\n",st7-st0);
   }

    /* -----------  Clean up ----------------- */   
   fftw_destroy_plan(p_before);
   fftw_destroy_plan(p_after);
   fftw_destroy_plan(p_veryend);

   ltfat_free(ff);
   ltfat_free(cf);
   
}
