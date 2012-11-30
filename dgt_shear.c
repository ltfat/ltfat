#include <complex.h>
#include "config.h"
#include <math.h>
#include "ltfat.h"
#include <stdio.h>

#define PI 3.1415926535897932384626433832795

LTFAT_EXTERN void
LTFAT_NAME(pchirp)(const int L, const int n, LTFAT_COMPLEX *g)
{ 

   const LTFAT_REAL LL=2.0*L;
   const LTFAT_REAL Lpone=L+1;
   
   for (int m=0;m<L;m++)
   {
      g[m] = cexp(I*PI*fmod(Lpone*n*m*m,LL)/L);
   }
  
}


LTFAT_EXTERN LTFAT_NAME(dgt_shear_plan)
LTFAT_NAME(dgt_shear_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
				 const int L, const int W, const int a, const int M,
				 const int s0, const int s1, const int br,
				 LTFAT_COMPLEX *cout,
				 unsigned flags)
{ 
   LTFAT_NAME(dgt_shear_plan) plan;
   
   plan.a=a;
   plan.M=M;
   plan.L=L;
   plan.W=W;

   plan.s0=s0;
   plan.s1=s1;
   plan.br=br;

   const int b=L/M;
   const int N=L/a;

   const int ar = a*b/br;
   const int Mr = L/br;
   const int Nr = L/ar;

   plan.f     = (LTFAT_COMPLEX *)f;
   plan.fwork = (LTFAT_COMPLEX *)f;
   plan.gwork = (LTFAT_COMPLEX *)g;
   plan.cout  = cout;

   plan.c_rect = ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEX));

   LTFAT_COMPLEX *f_before_fft = (LTFAT_COMPLEX *)f;
   LTFAT_COMPLEX *g_before_fft = (LTFAT_COMPLEX *)g;

   if ((!s0==0) || (!s1==0))
   {
      plan.fwork = ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
      plan.gwork = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));	 
   }

            
   if (!s1==0)
   {
      plan.p1 = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

      plan.fwork = ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
      plan.gwork = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

      LTFAT_NAME(pchirp)(L,s1,plan.p1);

      for (int l=0;l<L;l++)
      {
	 plan.gwork[l] = g[l]*plan.p1[l];
      }
      
      f_before_fft=plan.fwork;
      g_before_fft=plan.gwork;

   }
   
   if (!s0==0)
   {
      /* Allocate memory and compute the pchirp */
      plan.p0 = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));
      LTFAT_NAME(pchirp)(L,-s0,plan.p0);

      /* if data has already been copied to the working arrays, use
       * inline FFTs. Otherwise, if this is the first time they are
       * being used, do the copying using the fft. */

      plan.f_plan = LTFAT_FFTW(plan_many_dft)(1, &L, W,
					      f_before_fft, NULL, 1, L,
					      plan.fwork, NULL, 1, L,
					      FFTW_FORWARD, flags);

      plan.g_plan = LTFAT_FFTW(plan_dft_1d)(L, g_before_fft, plan.gwork, FFTW_FORWARD, flags);

      /* Execute the FFTs */
      LTFAT_FFTW(execute)(plan.g_plan);

      /* Multiply g by the chirp and scale by 1/L */
      for (int l=0;l<L;l++)
      {
	 plan.gwork[l] = plan.gwork[l]*plan.p0[l]/L;
      }
      
      /* Call the rectangular computation in the frequency domain*/
      /* LTFAT_NAME(dgt_long)(plan.fwork,plan.gwork,L,W,br,Nr,plan.c_rect); */
      /* Call the rectangular computation in the frequency domain*/
      plan.rect_plan = LTFAT_NAME(dgt_long_init)(plan.fwork, plan.gwork, L, W, br, Nr, plan.c_rect, flags);

   }
   else 
   {     

      /* Call the rectangular computation in the time domain */
      /* LTFAT_NAME(dgt_long)(plan.fwork,plan.gwork,L,W,ar,Mr,plan.c_rect); */

      plan.rect_plan = LTFAT_NAME(dgt_long_init)(plan.fwork, plan.gwork, 
						 L, W, ar, Mr, plan.c_rect, flags);      

   }  

   plan.finalmod = ltfat_malloc(2*N*sizeof(LTFAT_COMPLEX));
   
   for (int n=0;n<2*N;n++)
   {
      plan.finalmod[n]=cexp(PI*I*n/N);
   }

   return plan;

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shear_execute)(const LTFAT_NAME(dgt_shear_plan) plan)
{

   const int a=plan.a;
   const int M=plan.M;
   const int L=plan.L;

   const int b=plan.L/plan.M;
   const int N=plan.L/plan.a;
   const int s0 = plan.s0;
   const int s1 = plan.s1;

   const int ar = plan.a*b/plan.br;
   const int Mr = plan.L/plan.br;
   const int Nr = plan.L/ar;


   if (!plan.s1==0)
   {
      for (int w=0;w<plan.W;w++)
      {
      	 for (int l=0;l<plan.L;l++)
      	 {
	    plan.fwork[l+w*plan.L] = plan.f[l+w*plan.L]*plan.p1[l];
      	 }
      }
      
   }

   const int cc1=ar/a;
   const int cc2=-plan.s0*plan.br/a;
   const int twoN=2*N;
   const int cc6=(s0*s1+1)*plan.br;

   if (!plan.s0==0)
   {

      const int cc3=a*plan.s1*(L+1);
      const int cc4=cc2*plan.br*(L+1);
      const int cc5=2*cc1*plan.br;

      LTFAT_FFTW(execute)(plan.f_plan);

      for (int w=0;w<plan.W;w++)
      {
	 for (int l=0;l<plan.L;l++)
	 {

	    plan.fwork[l+w*plan.L] = plan.fwork[l+w*plan.L]*plan.p0[l];
	 }
      }

      LTFAT_NAME(dgt_long_execute)(plan.rect_plan);

      for (int k=0;k<Nr;k++)
      {
      	 for (int m=0;m<Mr;m++)
      	 {      	    
	    const int sq1 = k*cc1+cc2*m;
	    
	    const int phsidx = positiverem(cc3*sq1*sq1-m*(cc4*m+k*cc5),twoN);
                
	    const int idx1 = positiverem(k*cc1+cc2*m,N);
	    
	    /* The line below has a hidden floor operation by diving with the last b */
	    const int idx2 = positiverem(-s1*k*ar+cc6*m,L)/b;

	    for (int w=0;w<plan.W;w++) 
	    {
	       const int inidx  = positiverem(-k,Nr)+m*Nr+w*plan.M*N;
	       const int outidx = idx2+idx1*M+w*M*N;
	       plan.cout[outidx] = plan.c_rect[inidx]*plan.finalmod[phsidx];

	    }
      	 }
      }

   }
   else
   {
      const int cc3 = s1*(L+1);
      const int cc4 = s0*plan.br*plan.br*(L+1);

      LTFAT_NAME(dgt_long_execute)(plan.rect_plan);

      for (int k=0;k<Nr;k++)
      {
      	 for (int m=0;m<Mr;m++)
      	 {
	    const int sq1=k*ar-s0*m*plan.br;
	    const int phsidx= positiverem((cc3*sq1*sq1+cc4*m*m)/a,twoN);
                
	    const int idx1 = positiverem(k*cc1+cc2*m,N);
	    
	    /* The line below has a hidden floor operation by diving with the last b */
	    const int idx2 = positiverem(-s1*k*ar+cc6*m,L)/b;
            	    
      	    for (int w=0;w<plan.W;w++)
      	    {
      	       const int inidx  = m+k*Mr+w*M*N;
      	       const int outidx = idx2+idx1*M+w*M*N;
      	       plan.cout[outidx] = plan.c_rect[inidx]*plan.finalmod[phsidx];
      	    }
      	 }
      }

   }      
}
   

LTFAT_EXTERN void
LTFAT_NAME(dgt_shear_done)(LTFAT_NAME(dgt_shear_plan) plan)
{
   ltfat_free(plan.finalmod);

   LTFAT_NAME(dgt_long_done)(plan.rect_plan);
  
   ltfat_free(plan.c_rect);

   if (!plan.s0==0)
   {
      ltfat_free(plan.p0);
   }

   if (!plan.s1==0)
   {
      ltfat_free(plan.p1);
   }

   if ((!plan.s0==0) || (!plan.s1==0))
   {
      ltfat_free(plan.fwork);
      ltfat_free(plan.gwork);
   }
   
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_shear)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
		      const int L, const int W, const int a, const int M,
		      const int s0, const int s1, const int br,
		      LTFAT_COMPLEX *cout)
{ 

   LTFAT_NAME(dgt_shear_plan) plan = LTFAT_NAME(dgt_shear_init)(
      f,g,L,W,a,M,s0,s1,br,cout,FFTW_ESTIMATE);

   LTFAT_NAME(dgt_shear_execute)(plan);

   LTFAT_NAME(dgt_shear_done)(plan);

}
   

   

