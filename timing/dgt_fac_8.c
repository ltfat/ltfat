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
    conserve memory and hopefully use the cache hierachy better 

    The routine uses a very small buffer to do the DFTs.

    Integer indexing is optimized.

    Special code for integer oversampling.

    Special code for p=1 q=2, p=2, q=3
*/

void dgt_fac(ltfat_complex *f, ltfat_complex *gf, const int L, const int W,
	     const int R, const int a, const int M, ltfat_complex *cout, int dotime)
{

   /*  --------- initial declarations -------------- */

   int b, N, c, d, p, q, h_a, h_m;
   
   ltfat_complex *gbase, *fbase, *cbase;

   int l, k, r, s, u, w, rw, nm, mm, km;
   int ld2a, ld1b, ld3b, find, cfind;
   int ld4c, ld5c;
   int rem;

   fftw_plan p_before, p_after, p_veryend;
   ltfat_complex *ff, *cf, *sbuf, *fp, *ffp;

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
   sbuf = (ltfat_complex*)ltfat_malloc(d*sizeof(ltfat_complex));

   /* Create plans. In-place. */

   p_before = fftw_plan_dft_1d(d, sbuf, sbuf, FFTW_FORWARD,
				   FFTW_MEASURE);

   p_after  = fftw_plan_dft_1d(d, sbuf, sbuf, FFTW_BACKWARD,
				   FFTW_MEASURE);
   
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
      printf("DGT_FAC_7: Planning phase %f\n",st1-st0);
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
      find=0;
      ffp=ff;
      fp=f+r;
      if (p==1)

	 /* Integer oversampling case */
      {
	 
	 for (w=0;w<W;w++)
	 {
	    for (l=0;l<q;l++)
	    {
	       for (s=0;s<d;s++)
	       {		  
		  rem = positiverem(s*M+l*a, L);
		  
#ifdef HAVE_COMPLEX_H
		  sbuf[s]=fp[rem];
#else
		  sbuf[s][0]=fp[rem][0];
		  sbuf[s][1]=fp[rem][1];
#endif
		  
	       }
	       
	       fftw_execute(p_before);
	       
	       for (s=0;s<d;s++)
	       {		  
#ifdef HAVE_COMPLEX_H
		 ffp[s*ld2a]=sbuf[s]*scalconst;
#else
		 ffp[s*ld2a][0]=sbuf[s][0]*scalconst;
		 ffp[s*ld2a][1]=sbuf[s][1]*scalconst;
#endif		     
		  
	       }
	       ffp++;
	    }
	    fp+=L;
	 }
	 fp-=L*W;
	 
      }
      else
      {      
	 /* rational sampling case */

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
		     sbuf[s]=fp[rem];
#else
		     sbuf[s][0]=fp[rem][0];
		     sbuf[s][1]=fp[rem][1];
#endif
		     
		  }
		  
		  fftw_execute(p_before);
		  
		  for (s=0;s<d;s++)
		  {		  
#ifdef HAVE_COMPLEX_H
		     ffp[s*ld2a]=sbuf[s]*scalconst;
#else
		     ffp[s*ld2a][0]=sbuf[s][0]*scalconst;
		     ffp[s*ld2a][1]=sbuf[s][1]*scalconst;
#endif		     
		     
		  }
		  ffp++;
	       }
	    }
	    fp+=L;
	 }
	 fp-=L*W;
      }

      /* ----------- compute matrix multiplication ----------- */

      /* Do the matmul  */
      if (p==1)
      {
	 /* Integer oversampling case */
	 

	 /* Rational oversampling case */
	 for (s=0;s<d;s++)
	 {	
	    gbase=gf+(r+s*c)*q*R;
	    fbase=ff+s*q*W;
	    cbase=cf+s*q*q*W*R;
	    
	    for (nm=0;nm<q*W;nm++)
	    {
	       for (mm=0;mm<q*R;mm++)
	       {
#ifdef HAVE_COMPLEX_H
		  *cbase=conj((*gbase))*(*fbase);
#else
		  (*cbase)[0]=(*gbase)[0]*(*fbase)[0]+(*gbase)[1]*(*fbase)[1];
		  (*cbase)[1]=(*gbase)[0]*(*fbase)[1]-(*gbase)[1]*(*fbase)[0];
#endif
		  gbase++;
		  cbase++;
	       }			       
	       gbase-=q*R;
	       fbase++;
	    }
	    cbase-=q*R*q*W;
	 }




      }
      else
      {

	 /* Rational oversampling case */
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
		  *cbase=0.0;
		  for (km=0;km<p;km++)
		  {
		     *cbase+=conj((*gbase))*fbase[find+nm*p];
		     gbase+=1;
		     find+=1;
		  }
#else
		  (*cbase)[0]=0.0;
		  (*cbase)[1]=0.0;
		  for (km=0;km<p;km++)
		  {
		     (*cbase)[0]+=(*gbase)[0]*(*fbase)[0]+(*gbase)[1]*(*fbase)[1];
		     (*cbase)[1]+=(*gbase)[0]*(*fbase)[1]-(*gbase)[1]*(*fbase)[0];
		     gbase++;
		     fbase++;
		  }
#endif
		  fbase-=p;
		  cbase++;
	       }			       
	       gbase-=q*R*p;
	       fbase+=p;
	    }
	    cbase-=q*R*q*W;
	    fbase-=p*q*W;
	 }
      }



      /*  -------  compute inverse coefficient factorization ------- */
      cfind=0;
      ld4c=M*N;
      ld5c=M*N*R;

      /* Cover both integer and rational sampling case */
      for (w=0;w<W;w++)
      {
	 /* Complete inverse fac of coefficients */
	 for (l=0;l<q;l++)
	 {
	    for (rw=0;rw<R;rw++)
	    {
	       for (u=0;u<q;u++)
	       {	       	       
		  for (s=0;s<d;s++)	       
		  {	
#ifdef HAVE_COMPLEX_H	  
		     sbuf[s]   =cf[cfind+s*ld3b];
#else
		     sbuf[s][0]=cf[cfind+s*ld3b][0];
		     sbuf[s][1]=cf[cfind+s*ld3b][1];
#endif
		  }
		  cfind+=1;
		  
		  /* Do inverse fft of length d */
		  fftw_execute(p_after);
		  
		  for (s=0;s<d;s++)	       
		  {	
		     rem= r+l*c+positiverem(u+s*q-l*h_a,N)*M+rw*ld4c+w*ld5c;
#ifdef HAVE_COMPLEX_H	  
		     cout[rem]=sbuf[s];
#else
		     cout[rem][0]=sbuf[s][0];
		     cout[rem][1]=sbuf[s][1];
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
      printf("DGT_FAC_7: Main loop done %f\n",st6-st1);
   }

   /* FFT to modulate the coefficients. */
   fftw_execute(p_veryend);   

   if (dotime)
   {
      st7=ltfat_time();
      printf("DGT_FAC_7: Final FFT %f\n",st7-st6);
      printf("DGT_FAC_7: Total time %f\n",st7-st0);
   }

    /* -----------  Clean up ----------------- */   
   fftw_destroy_plan(p_before);
   fftw_destroy_plan(p_after);
   fftw_destroy_plan(p_veryend);

   ltfat_free(sbuf);
   ltfat_free(ff);
   ltfat_free(cf);
   
}

