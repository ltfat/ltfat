#include <complex.h>
#include "config.h"
#include <math.h>
#include "ltfat.h"
#include <stdio.h>

#define PI 3.1415926535897932384626433832795

static inline long positiverem_long(long int a,int b)
{
  const long c = a%b;
  return(c<0 ? c+b : c);
}


LTFAT_EXTERN void
LTFAT_NAME(pchirp)(const long L, const long n, LTFAT_COMPLEXH *g)
{

   const long LL=2*L;
   const long Lponen=positiverem_long((L+1)*n,LL);

   for (long m=0;m<L;m++)
   {
      const long idx = positiverem_long(
   	 positiverem_long(Lponen*m,LL)*m,LL);

      g[m] = cexp(1.0*I*PI*idx/L);
   }


   /* const LTFAT_REAL LL=2.0*L; */
   /* const LTFAT_REAL Lpone=L+1; */

   /* for (int m=0;m<L;m++) */
   /* { */
   /*    //g[m] = cexp(I*PI*fmod(Lpone*n*m*m,LL)/L); */
   /*    g[m] = cexp(I*PI*fmod(fmod(fmod(Lpone*n,LL)*m,LL)*m,LL)/L); */
   /* } */

}


LTFAT_EXTERN LTFAT_NAME(dgt_shear_plan)
LTFAT_NAME(dgt_shear_init)(const LTFAT_COMPLEXH *f, const LTFAT_COMPLEXH *g,
			   const int L, const int W, const int a, const int M,
			   const int s0, const int s1, const int br,
			   LTFAT_COMPLEXH *cout,
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

   plan.f     = (LTFAT_COMPLEXH *)f;
   plan.fwork = (LTFAT_COMPLEXH *)f;
   plan.gwork = (LTFAT_COMPLEXH *)g;
   plan.cout  = cout;

   plan.c_rect = ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEXH));

   LTFAT_COMPLEXH *f_before_fft = (LTFAT_COMPLEXH *)f;
   LTFAT_COMPLEXH *g_before_fft = (LTFAT_COMPLEXH *)g;

   if ((s0!=0) || (s1!=0))
   {
      plan.fwork = ltfat_malloc(L*W*sizeof(LTFAT_COMPLEXH));
      plan.gwork = ltfat_malloc(L*sizeof(LTFAT_COMPLEXH));
   }


   if (s1)
   {
      plan.p1 = ltfat_malloc(L*sizeof(LTFAT_COMPLEXH));

      LTFAT_NAME(pchirp)(L,s1,plan.p1);

      for (int l=0;l<L;l++)
      {
	 plan.gwork[l] = g[l]*plan.p1[l];
      }

      f_before_fft=plan.fwork;
      g_before_fft=plan.gwork;

   }

   if (s0==0)
   {

      /* Call the rectangular computation in the time domain */
      /* LTFAT_NAME(dgt_long)(plan.fwork,plan.gwork,L,W,ar,Mr,plan.c_rect); */

      plan.rect_plan = LTFAT_NAME(dgt_long_init)(plan.fwork, plan.gwork,
						 L, W, ar, Mr, plan.c_rect, flags);
   }
   else
   {

      /* Allocate memory and compute the pchirp */
      plan.p0 = ltfat_malloc(L*sizeof(LTFAT_COMPLEXH));
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

   plan.finalmod = ltfat_malloc(2*N*sizeof(LTFAT_COMPLEXH));

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
   const long s0 = plan.s0;
   const long s1 = plan.s1;

   const int ar = plan.a*b/plan.br;
   const int Mr = plan.L/plan.br;
   const int Nr = plan.L/ar;


   if (s1)
   {
      for (int w=0;w<plan.W;w++)
      {
      	 for (int l=0;l<plan.L;l++)
      	 {
	    plan.fwork[l+w*plan.L] = plan.f[l+w*plan.L]*plan.p1[l];
      	 }
      }

   }


   if (s0==0)
   {

      const int twoN=2*N;

      /* In this case, cc1=1 */

      const long cc3 = positiverem_long(s1*(L+1),twoN);

      const long tmp1=positiverem_long(cc3*a,twoN);

      LTFAT_NAME(dgt_long_execute)(plan.rect_plan);

      for (int k=0;k<N;k++)
      {
	 const int phsidx= positiverem_long((tmp1*k)%twoN*k,twoN);
	 const long part1= positiverem_long(-s1*k*a,L);
      	 for (int m=0;m<M;m++)
      	 {
	    /* The line below has a hidden floor operation when dividing with the last b */
	    const int idx2 = ((part1+b*m)%L)/b;

	    const int inidx  =    m+k*M;
	    const int outidx = idx2+k*M;
      	    for (int w=0;w<plan.W;w++)
      	    {
      	       plan.cout[outidx+w*M*N] = plan.c_rect[inidx+w*M*N]*plan.finalmod[phsidx];
      	    }
      	 }
      }


   }
   else
   {

      const int twoN=2*N;
      const long cc1=ar/a;
      const long cc2=positiverem_long(-s0*plan.br/a,twoN);
      const long cc3=positiverem_long(a*s1*(L+1),twoN);
      const long cc4=positiverem_long(cc2*plan.br*(L+1),twoN);
      const long cc5=positiverem_long(2*cc1*plan.br,twoN);
      const long cc6=positiverem_long((s0*s1+1)*plan.br,L);

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
	 const long part1 = positiverem_long(-s1*k*ar,L);
      	 for (int m=0;m<Mr;m++)
      	 {
	    const long sq1 = k*cc1+cc2*m;

	    const int phsidx = positiverem_long(
	       (cc3*sq1*sq1)%twoN-(m*(cc4*m+k*cc5))%twoN,twoN);

	    /* The line below has a hidden floor operation when dividing with the last b */
	    const int idx2 = ((part1+cc6*m)%L)/b;

	    const int inidx  = positiverem(-k,Nr)+m*Nr;
	    const int outidx = idx2+(sq1%N)*M;
	    for (int w=0;w<plan.W;w++)
	    {
	       plan.cout[outidx+w*M*N] = plan.c_rect[inidx+w*M*N]*plan.finalmod[phsidx];

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

   if (plan.s0 || plan.s1)
   {
      ltfat_free(plan.fwork);
      ltfat_free(plan.gwork);
   }

   if (plan.s0)
   {
      ltfat_free(plan.p0);
   }

   if (plan.s1)
   {
      ltfat_free(plan.p1);
   }


}


LTFAT_EXTERN void
LTFAT_NAME(dgt_shear)(const LTFAT_COMPLEXH *f, const LTFAT_COMPLEXH *g,
		      const int L, const int W, const int a, const int M,
		      const int s0, const int s1, const int br,
		      LTFAT_COMPLEXH *cout)
{

   LTFAT_NAME(dgt_shear_plan) plan = LTFAT_NAME(dgt_shear_init)(
      f,g,L,W,a,M,s0,s1,br,cout,FFTW_ESTIMATE);

   LTFAT_NAME(dgt_shear_execute)(plan);

   LTFAT_NAME(dgt_shear_done)(plan);

}




