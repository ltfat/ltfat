#include "config.h"
#include <math.h>
#include "ltfat.h"

#define PI 3.1415926535897932384626433832795

LTFAT_EXTERN void
LTFAT_NAME(pchirp)(const int L, const int n, LTFAT_COMPLEX *g)
{ 

   const double LL=2.0*L;
   const double Lpone=L+1;
   
   for (int m=0;m<L;m++)
   {
      const double work = PI*fmod(Lpone*n*m*m,LL)/L;
      g[m][0] = cos(work);
      g[m][1] = sin(work);
   }
  
}


LTFAT_EXTERN void
LTFAT_NAME(nonsepwin2multi)(const LTFAT_COMPLEX *g,
			    const int L, const int a, const int M,
			    const int lt1, const int lt2,
			    LTFAT_COMPLEX *mwin)
{ 

  const int b=L/M;

  const LTFAT_REAL scal = 2*PI/L;

  for (int w=0;w<lt2;w++)
  {
     const int wavenum=((w*lt1)%lt2)*b/lt2;
     for (int l=0;l<L;l++)
     {
	const LTFAT_REAL e0 = cos(scal*l*wavenum);
	const LTFAT_REAL e1 = sin(scal*l*wavenum);
	const int idx = positiverem(l-w*a,L);
	mwin[l+w*L][0]=e0*g[idx][0]-e1*g[idx][1];
	mwin[l+w*L][1]=e1*g[idx][0]+e0*g[idx][1];
     }
  }
}


LTFAT_EXTERN void
LTFAT_NAME(nonsepdgt_multi)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			    const int L, const int W, const int a, const int M,
			    const int lt1, const int lt2,
			    LTFAT_COMPLEX *c)
{ 
   const int N   = L/a;
   const int Ns  = N/lt2;
   LTFAT_COMPLEX *mwin = ltfat_malloc(L*lt2*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *c_scratch = ltfat_malloc(M*Ns*W*sizeof(LTFAT_COMPLEX));

   LTFAT_NAME(nonsepwin2multi)(g,L,a,M,lt1,lt2,mwin);
   
   for (int win=0;win<lt2;win++)
   {
      LTFAT_NAME(dgt_long)(f,mwin+L*win,L,W,a*lt2,M,c_scratch);
      for (int n=0;n<Ns;n++)
      {
   	 const LTFAT_REAL scal  = -2*PI*(a*n*((win*lt1)%lt2))/M;
   	 const LTFAT_REAL scal0 = cos(scal);
   	 const LTFAT_REAL scal1 = sin(scal);
   	 for (int w=0;w<W;w++)
   	 {
   	    for (int m=0;m<M;m++)
   	    {
	       const int idx_s = m+n*M+w*M*Ns;
	       const int idx = m+win*M+n*lt2*M+w*M*N;
	       c[idx][0] = scal0*c_scratch[idx_s][0]-scal1*c_scratch[idx_s][1];
   	       c[idx][1] = scal1*c_scratch[idx_s][0]+scal0*c_scratch[idx_s][1];
   	    }
   	 }
      }
   }

   ltfat_free(c_scratch);
   ltfat_free(mwin);

}
