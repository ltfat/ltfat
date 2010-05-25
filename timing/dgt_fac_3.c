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


/*  This routine computes the DGT factorization using contiguous FFTs
    and reshuffles data in between the FFTs and the matrix product. */


void dgt_fac(ltfat_complex *f, ltfat_complex *gf, const int L, const int W,
	     const int R, const int a, const int M, ltfat_complex *cout, int dotime)
{

   /*  --------- initial declarations -------------- */

   int b, N, c, d, p, q, h_a, h_m;
   
   ltfat_complex *gbase, *fbase, *cbase;

   int l, k, r, s, u, w, rw, nm, mm, km;
   int ld2, ld3, ld2n, ld3n, ldc1, ldc2, ldc3, ldc1n, ldc2n, ldc3n, ldo2, ldo3;
   int rem;

   fftw_plan p_before, p_after, p_veryend;
   ltfat_complex  *ff, *cf, *tmp_ff, *tmp_cf;

   double scalconst;

   double st0, st1, st2, st3, st4, st5, st6, st7, st8, st9;
   
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

   ff     = (ltfat_complex*)ltfat_malloc(L*W*sizeof(ltfat_complex));
   tmp_ff = (ltfat_complex*)ltfat_malloc(L*W*sizeof(ltfat_complex));
   cf     = (ltfat_complex*)ltfat_malloc(c*d*q*q*W*R*sizeof(ltfat_complex));
   tmp_cf = (ltfat_complex*)ltfat_malloc(c*d*q*q*W*R*sizeof(ltfat_complex));

   /* Create plans. In-place. Contiguous. */
   
   p_before = fftw_plan_many_dft(1, &d, c*p*q*W,
				 tmp_ff, NULL,
				 1, d,
				 tmp_ff, NULL,
				 1, d,
				 FFTW_FORWARD, FFTW_OPTITYPE);


   p_after = fftw_plan_many_dft(1, &d, c*q*q*W*R,
				tmp_cf, NULL,
				1, d,
				tmp_cf, NULL,
				1, d,
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
      printf("DGT_FAC_3: Planning phase %f\n",st1-st0);
   }   

   /*  ---------- compute signal factorization ----------- */

   /* Leading dimensions of the 4dim array 'ff'. */
   ld2=p*q*W;
   ld3=c*p*q*W;

   /* Leading dimensions of the temporary array 'tmp_ff' suited for contiguous
    * FFTs.
    */
   ld2n=d*p;
   ld3n=d*p*q*W;

   /* Leading dimensions of cf */
   ldc1=q*R;
   ldc2=q*R*q*W;
   ldc3=c*q*R*q*W;

   /* Leading dimensions of the temporary array 'tmp_cf' suited for contiguous
    * FFTs.
    */
   ldc1n=d*q;
   ldc2n=d*q*R;
   ldc3n=d*q*q*W*R;

   /* Leading dimensions of cout */
   ldo2=M*N;
   ldo3=M*N*R;

   for (r=0;r<c;r++)
   {
      for (w=0;w<W;w++)
      {
	 for (l=0;l<q;l++)
	 {
	   for (k=0;k<p;k++)
	   {	   
	      for (s=0;s<d;s++)
	      {
		 rem = positiverem(k*M+s*p*M-l*h_a*a, L);	 
#ifdef HAVE_COMPLEX_H
		 tmp_ff[s+k*d+(l+q*w)*ld2n+r*ld3n]   =f[r+rem+L*w]*scalconst;
#else
		 tmp_ff[s+k*d+(l+q*w)*ld2n+r*ld3n][0]=f[r+rem+L*w][0]*scalconst;
		 tmp_ff[s+k*d+(l+q*w)*ld2n+r*ld3n][1]=f[r+rem+L*w][1]*scalconst;
#endif
	      }
	    }
	 }
      }
   }           

   if (dotime)
   {
      st2=ltfat_time();
      printf("DGT_FAC_3: First permutation %f\n",st2-st1);
   }

   
   /* fft */
   fftw_execute(p_before);

   if (dotime)
   {
      st3=ltfat_time();
      printf("DGT_FAC_3: First FFT %f\n",st3-st2);
   }

   /* Do dimension permutation to complete signal factorization.
    * The k, r, w and l loop have been condensed into the same loop.
    */
   for (s=0;s<d;s++)
   {
     for (k=0;k<p*q*W*c;k++)
     {
#ifdef HAVE_COMPLEX_H
	ff[k+s*ld3]   =tmp_ff[s+k*d];
#else
	ff[k+s*ld3][0]=tmp_ff[s+k*d][0];
	ff[k+s*ld3][1]=tmp_ff[s+k*d][1];
#endif
      }
   }

   if (dotime)
   {
      st4=ltfat_time();
      printf("DGT_FAC_3: First reshuffle %f\n",st4-st3);
   }
   

   /* ----------- compute matrix multiplication ----------- */


   /* Do the matmul  */
   for (r=0;r<c;r++)
   {
      for (s=0;s<d;s++)
      {	
	 gbase=gf+(r+s*c)*p*q*R;
	 fbase=ff+(r+s*c)*p*q*W;
	 cbase=cf+(r+s*c)*q*q*W*R;

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
	       /* Scale because of FFTWs normalization. */
#endif
	    }		  
	 }	      	 
      }
   }

   if (dotime)
   {
      st5=ltfat_time();
      printf("DGT_FAC_3: Matrix multiplication %f\n",st5-st4);
   }


   /*  -------  compute inverse coefficient factorization ------- */

   /* Do the dimension permutation of cf. 
    * The r, w and l loop have been condensed into the same loop.
    */
   for (k=0;k<q*q*R*W*c;k++)
   {
      for (s=0;s<d;s++)
      {
#ifdef HAVE_COMPLEX_H
	 tmp_cf[s+k*d]   =cf[k+s*ldc3];
#else
	 tmp_cf[s+k*d][0]=cf[k+s*ldc3][0];
	 tmp_cf[s+k*d][1]=cf[k+s*ldc3][1];
#endif
      }
   }

   if (dotime)
   {
      st6=ltfat_time();
      printf("DGT_FAC_3: Second reshuffle %f\n",st6-st4);
   }

   /* Do inverse fft of length d */
   fftw_execute(p_after);
         
   if (dotime)
   {
      st7=ltfat_time();
      printf("DGT_FAC_3: Second FFT: IFFT %f\n",st7-st6);
   }

   /* Complete inverse fac of coefficients */
   for (rw=0;rw<R;rw++)
   {
      for (w=0;w<W;w++)
      {
	 for (s=0;s<d;s++)
	 {
	    for (u=0;u<q;u++)
	    {	       
	       for (l=0;l<q;l++)
	       {
		  rem = positiverem(u+s*q-l*h_a,N)*M;
		  for (r=0;r<c;r++)
		  {	
#ifdef HAVE_COMPLEX_H	  
		     cout[r+l*c+rem+rw*ldo2+w*ldo3]   =tmp_cf[s+u*d+rw*ldc1n+(l+q*w)*ldc2n+r*ldc3n];
#else
		     cout[r+l*c+rem+rw*ldo2+w*ldo3][0]=tmp_cf[s+u*d+rw*ldc1n+(l+q*w)*ldc2n+r*ldc3n][0];
		     cout[r+l*c+rem+rw*ldo2+w*ldo3][1]=tmp_cf[s+u*d+rw*ldc1n+(l+q*w)*ldc2n+r*ldc3n][1];
#endif
		  }
	       }
	    }
	 }
      }      
   }     

   if (dotime)
   {
      st8=ltfat_time();
      printf("DGT_FAC_3: Second permutation %f\n",st8-st7);
   }

   /* FFT to modulate the coefficients. */
   fftw_execute(p_veryend);   

   if (dotime)
   {
      st9=ltfat_time();
      printf("DGT_FAC_3: Final FFT %f\n",st9-st8);
      printf("DGT_FAC_3: Total time %f\n",st9-st0);
   }


    /* -----------  Clean up ----------------- */   

   fftw_destroy_plan(p_before);
   fftw_destroy_plan(p_after);
   fftw_destroy_plan(p_veryend);

   ltfat_free(tmp_ff);
   ltfat_free(ff);
   ltfat_free(tmp_cf);
   ltfat_free(cf);
   
}
