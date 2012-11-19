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
      const LTFAT_REAL work = PI*fmod(Lpone*n*m*m,LL)/L;
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



LTFAT_EXTERN void
LTFAT_NAME(nonsepdgt_shear)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			    const int L, const int W, const int a, const int M,
			    const int s0, const int s1, const int br,
			    LTFAT_COMPLEX *c)
{ 
   const int b=L/M;
   const int N=L/a;

   const int ar = a*b/br;
   const int Mr = L/br;
   const int Nr = L/ar;

   int using_fwork = 0;
   int using_gwork = 0;

   const unsigned flags = FFTW_ESTIMATE;

   LTFAT_COMPLEX *fwork = (LTFAT_COMPLEX *)f;
   LTFAT_COMPLEX *gwork = (LTFAT_COMPLEX *)g;

            
   if (!s1==0)
   {
      LTFAT_COMPLEX *p = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

      using_fwork=1;
      using_gwork=1;

      fwork = ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
      gwork = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

      LTFAT_NAME(pchirp)(L,s1,p);

      for (int l=0;l<L;l++)
      {
	 gwork[l][0] = g[l][0]*p[l][0]-g[l][1]*p[l][1];
	 gwork[l][1] = g[l][1]*p[l][0]+g[l][0]*p[l][1];

      }

      for (int w=0;w<W;w++)
      {
	 for (int l=0;l<L;l++)
	 {
	    fwork[l+w*L][0] = f[l+w*L][0]*p[l][0]-f[l+w*L][1]*p[l][1];
	    fwork[l+w*L][1] = f[l+w*L][1]*p[l][0]+f[l+w*L][0]*p[l][1];	    
	 }
      }

      ltfat_free(p);
      
   }
   
   if (!s0==0)
   {
      /* Allocate memory and compute the pchirp */
      LTFAT_COMPLEX *p = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));
      LTFAT_NAME(pchirp)(L,-s0,p);


      LTFAT_FFTW(plan) f_plan, g_plan;

      /* if data has already been copied to the working arrays, use
       * inline FFTs. Otherwise, if this is the first time they are
       * being used, do the copying using the fft. */

      if (!using_fwork)
      { 
	 fwork = ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
	 
	 f_plan = LTFAT_FFTW(plan_many_dft)(1, &L, W,
					    f, NULL, 1, L,
					    fwork, NULL, 1, L,
					    FFTW_FORWARD, flags);
	 
	 using_fwork=1;      
      }
      else
      {
	 f_plan = LTFAT_FFTW(plan_many_dft)(1, &L, W,
					    fwork, NULL, 1, L,
					    fwork, NULL, 1, L,
					    FFTW_FORWARD, flags);
      }
      
      if (!using_gwork)
      {
	 gwork = ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

	 g_plan = LTFAT_FFTW(plan_dft_1d)(L, g, gwork, FFTW_FORWARD, flags);
	 
	 using_gwork=1;
      }
      else
      {
	 g_plan = LTFAT_FFTW(plan_dft_1d)(L, gwork, gwork, FFTW_FORWARD, flags);
      }


      /* Execute the FFTs */
      LTFAT_FFTW(execute)(f_plan);
      LTFAT_FFTW(execute)(g_plan);

      /* Multiply g by the chirp and scale by 1/L */
      for (int l=0;l<L;l++)
      {
	 const LTFAT_REAL tmp = (gwork[l][0]*p[l][0]-gwork[l][1]*p[l][1])/L; 	 
	 gwork[l][1] = (gwork[l][1]*p[l][0]+gwork[l][0]*p[l][1])/L;
	 gwork[l][0] = tmp;
      }

      for (int w=0;w<W;w++)
      {
	 for (int l=0;l<L;l++)
	 {
	    const LTFAT_REAL tmp = fwork[l+w*L][0]*p[l][0]-fwork[l+w*L][1]*p[l][1];
	    
	    fwork[l+w*L][1] = fwork[l+w*L][1]*p[l][0]+fwork[l+w*L][0]*p[l][1];	    
	    fwork[l+w*L][0] = tmp;
	 }
      }

      ltfat_free(p);

      LTFAT_COMPLEX *c_rect = ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEX));
      
      /* Call the rectangular computation in the frequency domain*/
      LTFAT_NAME(dgt_long)(fwork,gwork,L,W,br,Nr,c_rect);
      
      for (int k=0;k<Nr;k++)
      {   
	 for (int m=0;m<Mr;m++)
	 {
	    const int t1 = k*ar-s0*m*br;
	    const int t2 = m*br;

	    const LTFAT_REAL phs = 
	       PI*positiverem((s1*t1*t1+s0*t2*t2)*(L+1)-2*(k*ar*m*br),2*L)/L;
	    
	    const LTFAT_REAL phs0 = cos(phs);
	    const LTFAT_REAL phs1 = sin(phs);
            
	    const int idx1 =       positiverem(    k*ar       -s0*m*br,L)/a;
	    const int idx2 = floor(positiverem(-s1*k*ar+(s0*s1+1)*m*br,L)/b);
            
	    for (int w=0;w<W;w++)
	    {                  
	       const int inidx  = positiverem(-k,Nr)+m*Nr+w*M*N;
	       const int outidx = idx2+idx1*M+w*M*N;
	       c[outidx][0] = c_rect[inidx][0]*phs0-c_rect[inidx][1]*phs1;
	       c[outidx][1] = c_rect[inidx][1]*phs0+c_rect[inidx][0]*phs1;
	    }
	 }
      }

      ltfat_free(c_rect);
   }
   else 
   {     

      LTFAT_COMPLEX *c_rect = ltfat_malloc(M*N*W*sizeof(LTFAT_COMPLEX));
      
      /* Call the rectangular computation in the time domain */
      LTFAT_NAME(dgt_long)(fwork,gwork,L,W,ar,Mr,c_rect);

      for (int k=0;k<Nr;k++)
      {   
	 for (int m=0;m<Mr;m++)
	 {
	    const int t1 = k*ar-s0*m*br;
	    const int t2 = m*br;
	    const LTFAT_REAL phs = 
	       PI*positiverem((s1*t1*t1+s0*t2*t2)*(L+1),2*L)/L;
	    
	    const LTFAT_REAL phs0 = cos(phs);
	    const LTFAT_REAL phs1 = sin(phs);
            
	    const int idx1 =       positiverem(    k*ar       -s0*m*br,L)/a;
	    const int idx2 = floor(positiverem(-s1*k*ar+(s0*s1+1)*m*br,L)/b);
            
	    
	    for (int w=0;w<W;w++)
	    {
	       const int inidx  = m+k*Mr+w*M*N;
	       const int outidx = idx2+idx1*M+w*M*N;
	       c[outidx][0] = c_rect[inidx][0]*phs0-c_rect[inidx][1]*phs1;
	       c[outidx][1] = c_rect[inidx][1]*phs0+c_rect[inidx][0]*phs1;
	    }
	 }
      }

      ltfat_free(c_rect);
   }  

   if (using_fwork)
      ltfat_free(fwork);

   if (using_gwork)
      ltfat_free(gwork);

}



