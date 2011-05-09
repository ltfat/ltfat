#include <stdlib.h>
#include "config.h"
#include "fftw3.h"
#include "ltfat.h"

LTFAT_EXTERN
void LTFAT_NAME(dgt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			 const int L, const int W, const int a,
			 const int M, LTFAT_COMPLEX *cout)
{
  
  LTFAT_NAME(dgt_long_plan) plan =
    LTFAT_NAME(dgt_long_init)(f, g, L, W, a, M, cout, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgt_long_execute)(plan);

  LTFAT_NAME(dgt_long_done)(plan);
  
}

LTFAT_EXTERN
void LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			 const int L, const int W, const int a,
			 const int M, LTFAT_COMPLEX *cout)
{
 
  LTFAT_NAME(dgtreal_long_plan) plan =
    LTFAT_NAME(dgtreal_long_init)(f, g, L, W, a, M, cout, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgtreal_long_execute)(plan);

  LTFAT_NAME(dgtreal_long_done)(plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
		     const int L, const int gl,
		     const int W,  const int a, const int M, 
		     LTFAT_COMPLEX *cout)
{
 
  LTFAT_NAME(dgt_fb_plan) plan =
    LTFAT_NAME(dgt_fb_init)(g, gl, a, M, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgt_fb_execute)(plan, f, L, W, cout);

  LTFAT_NAME(dgt_fb_done)(plan);

}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
	      const int L, const int gl,
	      const int W, const int a, const int M, 
	      LTFAT_COMPLEX *cout)
{
 
  LTFAT_NAME(dgtreal_fb_plan) plan =
    LTFAT_NAME(dgtreal_fb_init)(g, gl, a, M, FFTW_ESTIMATE);
  
  LTFAT_NAME(dgtreal_fb_execute)(plan, f, L, W, cout);

  LTFAT_NAME(dgtreal_fb_done)(plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			const int L, const int gl,
			const int W, const int a, const int M, const int bl, 
			LTFAT_COMPLEX *cout)
{

   const int N       = L/a;
   const int M2      = M/2+1;
   const int Lext    = bl+gl;
   const int Nb      = L/bl;
   const int b2      = gl/a/2;
   const int Nblock  = bl/a;
   const int Nblocke = Lext/a;

   const LTFAT_REAL *buf = ltfat_malloc(Lext*W*sizeof(LTFAT_REAL));

   const LTFAT_REAL *gext = ltfat_malloc(Lext*sizeof(LTFAT_REAL));

 
   const LTFAT_COMPLEX *cbuf = ltfat_malloc(M2*W*Lext/a*sizeof(LTFAT_COMPLEX));
 
   LTFAT_NAME(fir2iir_r)(g, gl, Lext, gext);

   /* Zero the output array, as we will be adding to it */
   for (ii=0; ii<M*N*W; ii++}
   {
      cout[ii][0]=0.0;
      cout[ii][1]=0.0;
   }

   LTFAT_NAME(dgtreal_long_plan) plan =
      LTFAT_NAME(dgtreal_long_init)(buf, gext, bl+gl, W, a, M,
				    cbuf, FFTW_ESTIMATE);
   
   for (int ii=0; ii<Nb; i++)
   {
      int s_ii;

      /* Copy to working buffer and zero the last part. */
      for (int w=0; w<W; w++)
      {
	 memcpy(f+w*L,buf+Lext*w,sizeof(LTFAT_REAL)*bl);
	 for (int jj=bl; jj<Lext;j++)
	 {
	    f[jj]=0.0;
	 }
      }
      

      for (int w=0; w<W; w++)
      {
	 LTFAT_COMPLEX coef_p;
	 LTFAT_COMPLEX cbuf_p;

	 /* Place large block */
	 coef_p = coef + w*M*N ;
	 cbuf_p = cbuf + w*M*Nblocke;

	 for (int m=0; m<M; m++)
	 {
	    for (int n=0;n<Nblock;n++)
	    {
	       coef_p[m+n*M][0] += cbuf_p[m+n*M][0];
	       coef_p[m+n*M][1] += cbuf_p[m+n*M][1];
	    }
	 }

	 /* Small block + */
	 s_ii=positiverem(ii+1,Nb);
	 coef(:,s_ii*Nblock+1   :s_ii*Nblock+b2)+=block(:,Nblock+1:Nblock+b2); 
	 for (int m=0; m<M; m++)
	 {
	    for (int n=0;n<Nblock;n++)
	    {
	       coef_p[m+n*M][0] += cbuf_p[m+n*M][0];
	       coef_p[m+n*M][1] += cbuf_p[m+n*M][1];
	    }
	 }

	 
	 /* Small block - */
	 s_ii=positiverem(ii-1,Nb);
	 coef(:,s_ii*Nblock-b2+1:s_ii*Nblock)   +=block(:,Nblock+b2+1:Nblock+2*b2);
	 
      }
       


      LTFAT_NAME(dgtreal_long_execute)(plan);
      
   }

   
   LTFAT_NAME(dgtreal_long_done)(plan);

   ltfat_free(cbuf);
   ltfat_free(gext);
   ltfat_free(buf);

}


