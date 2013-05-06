#include <complex.h>
#include "config.h"
#include <math.h>
#include "ltfat.h"
#include <stdio.h>

#define PI 3.1415926535897932384626433832795


LTFAT_EXTERN void
LTFAT_NAME(nonsepwin2multi)(const LTFAT_COMPLEXH *g,
			    const int L, const int Lg, const int a, const int M,
			    const int lt1, const int lt2,
			    LTFAT_COMPLEXH *mwin)
{

  const int b=L/M;

  const LTFAT_REAL scal = 2*PI/L;

  LTFAT_COMPLEXH *gwork = (LTFAT_COMPLEXH *)ltfat_malloc(L*sizeof(LTFAT_COMPLEXH));
  LTFAT_NAME(fir2long_c)((const LTFAT_COMPLEX*)g,Lg,L,(LTFAT_COMPLEX*)gwork);

  for (int w=0;w<lt2;w++)
  {
     const int wavenum=((w*lt1)%lt2)*b/lt2;
     for (int l=0;l<L;l++)
     {
 	mwin[l+w*L]=cexp(I*scal*l*wavenum)*gwork[positiverem(l-w*a,L)];
     }
  }

  ltfat_free(gwork);
}


LTFAT_EXTERN LTFAT_NAME(dgt_multi_plan)
LTFAT_NAME(dgt_multi_init)(const LTFAT_COMPLEXH *f, const LTFAT_COMPLEXH *g,
			   const int L, const int Lg, const int W, const int a, const int M,
			   const int lt1, const int lt2,
			   LTFAT_COMPLEXH *cout,unsigned flags)
{

   LTFAT_NAME(dgt_multi_plan) plan;

   plan.a=a;
   plan.M=M;
   plan.L=L;
   plan.Lg=Lg;
   plan.W=W;

   plan.lt1=lt1;
   plan.lt2=lt2;

   plan.f     = (LTFAT_COMPLEXH *)f;
   plan.cout  = cout;

   const int N   = L/a;
   const int Ns  = N/lt2;

   plan.mwin = (LTFAT_COMPLEXH *)ltfat_malloc(L*lt2*sizeof(LTFAT_COMPLEXH));
   plan.c_scratch = (LTFAT_COMPLEXH *)ltfat_malloc(M*Ns*W*sizeof(LTFAT_COMPLEXH));


   LTFAT_NAME(nonsepwin2multi)(g,L,Lg,a,M,lt1,lt2,plan.mwin);

   plan.rect_plan_array = (LTFAT_NAME(dgt_long_plan)*) ltfat_malloc(lt2*sizeof(LTFAT_NAME(dgt_long_plan)));

   for (int win=0;win<lt2;win++)
   {
      plan.rect_plan_array[win]=LTFAT_NAME(dgt_long_init)((const LTFAT_COMPLEX*)plan.f,(const LTFAT_COMPLEX*)(plan.mwin+L*win),L,W,a*lt2,M,
							 (LTFAT_COMPLEX*)plan.c_scratch,flags);
   }

   plan.mod = (LTFAT_COMPLEXH*) ltfat_malloc(N*sizeof(LTFAT_COMPLEXH));

   for (int win=0;win<plan.lt2;win++)
   {
      for (int n=0;n<Ns;n++)
      {
	   plan.mod[win+n*lt2] = cexp(-2*PI*I*(a*n*((win*lt1)%lt2))/M);
      }
   }

   return plan;
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_multi_execute)(const LTFAT_NAME(dgt_multi_plan) plan)
{
   const int N   = plan.L/plan.a;
   const int Ns  = N/plan.lt2;

   const int M=plan.M;
   const int W=plan.W;
   const int lt2=plan.lt2;

   for (int win=0;win<plan.lt2;win++)
   {
      LTFAT_NAME(dgt_long_execute)(plan.rect_plan_array[win]);
      for (int w=0;w<W;w++)
      {
	 for (int n=0;n<Ns;n++)
	 {
   	    for (int m=0;m<M;m++)
   	    {
	       plan.cout[m+win*M+n*lt2*M+w*M*N] = plan.mod[win+n*lt2]*plan.c_scratch[m+n*M+w*M*Ns];
   	    }
   	 }
      }
   }
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_multi_done)(LTFAT_NAME(dgt_multi_plan) plan)
{

   ltfat_free(plan.mod);

   for (int ii=0;ii<plan.lt2;ii++)
   {
      LTFAT_NAME(dgt_long_done)(plan.rect_plan_array[ii]);
   }

   ltfat_free(plan.rect_plan_array);

   ltfat_free(plan.c_scratch);
   ltfat_free(plan.mwin);
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_multi)(const LTFAT_COMPLEXH *f, const LTFAT_COMPLEXH *g,
		      const int L, const int Lg, const int W, const int a, const int M,
		      const int lt1, const int lt2,
		      LTFAT_COMPLEXH *cout)
{

   LTFAT_NAME(dgt_multi_plan) plan = LTFAT_NAME(dgt_multi_init)(
      f,g,L,Lg,W,a,M,lt1,lt2,cout,FFTW_ESTIMATE);

   LTFAT_NAME(dgt_multi_execute)(plan);

   LTFAT_NAME(dgt_multi_done)(plan);

}
