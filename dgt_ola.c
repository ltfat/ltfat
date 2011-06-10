#include <string.h>
#include "config.h"
#include "fftw3.h"
#include "ltfat.h"


LTFAT_EXTERN LTFAT_NAME(dgt_ola_plan)
LTFAT_NAME(dgt_ola_init)(const LTFAT_COMPLEX *g, const int gl,
			 const int W, const int a, const int M, const int bl,
			 unsigned flags)
{

   LTFAT_NAME(dgt_ola_plan) plan;

   plan.bl = bl;
   plan.gl = gl;
   plan.W  = W;

   const int Lext    = bl+gl;
   const int Nblocke = Lext/a;
   
   plan.buf  = ltfat_malloc(Lext*W*sizeof(LTFAT_COMPLEX));
   plan.gext = ltfat_malloc(Lext*sizeof(LTFAT_COMPLEX));   
   plan.cbuf = ltfat_malloc(M*Nblocke*W*sizeof(LTFAT_COMPLEX));
   
   LTFAT_NAME(fir2long_c)(g, gl, Lext, plan.gext);
   
   /* Zero the last part of the buffer, it will always be zero. */
   for (int w=0; w<W; w++)
   {      
      for (int jj=bl; jj<Lext;jj++)
      {
	 plan.buf[jj+w*Lext][0]=0.0;
	 plan.buf[jj+w*Lext][1]=0.0;
      }
   }
   
   plan.plan =
      LTFAT_NAME(dgt_long_init)((const LTFAT_COMPLEX*)plan.buf,
				(const LTFAT_COMPLEX*)plan.gext,
				Lext, W, a, M,
				plan.cbuf, flags);
   
   return (plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_execute)(const LTFAT_NAME(dgt_ola_plan) plan,
			    const LTFAT_COMPLEX *f, const int L,
			    LTFAT_COMPLEX *cout)

{
   const int bl      = plan.bl;
   const int gl      = plan.gl;
   const int a       = plan.plan.a;
   const int M       = plan.plan.M;
   const int N       = L/a;
   const int Lext    = bl+gl;
   const int Nb      = L/bl;
   const int b2      = gl/a/2;
   const int Nblock  = bl/a;
   const int Nblocke = Lext/a;
   const int W       = plan.W;


   /* Zero the output array, as we will be adding to it */
   for (int ii=0; ii<M*N*W; ii++)
   {
      cout[ii][0]=0.0;
      cout[ii][1]=0.0;
   }
   
   for (int ii=0; ii<Nb; ii++)
   {
      int s_ii;

      /* Copy to working buffer. */
      for (int w=0; w<W; w++)
      {
	 memcpy(plan.buf+Lext*w,f+ii*bl+w*L,sizeof(LTFAT_COMPLEX)*bl);
      }
      
      /* Execute the short DGT */
      LTFAT_NAME(dgt_long_execute)(plan.plan);

      /* Place the results */
      for (int w=0; w<W; w++)
      {
	 LTFAT_COMPLEX *cout_p;
	 LTFAT_COMPLEX *cbuf_p;

	 /* Place large block */
	 cout_p = cout + ii*M*Nblock+w*M*N ;
	 cbuf_p = plan.cbuf +             w*M*Nblocke;
	 for (int m=0; m<M; m++)
	 {
	    for (int n=0;n<Nblock;n++)
	    {
	       cout_p[m+n*M][0] += cbuf_p[m+n*M][0];
	       cout_p[m+n*M][1] += cbuf_p[m+n*M][1];
	    }
	 }

	 /* Small block + */
	 s_ii=positiverem(ii+1,Nb);
	 cout_p = cout + s_ii*M*Nblock+w*M*N ;
	 cbuf_p = plan.cbuf +      M*Nblock+w*M*Nblocke;
	 for (int m=0; m<M; m++)
	 {
	    for (int n=0;n<b2;n++)
	    {
	       cout_p[m+n*M][0] += cbuf_p[m+n*M][0];
	       cout_p[m+n*M][1] += cbuf_p[m+n*M][1];
	    }
	 }

	 
	 /* Small block - */
	 s_ii=positiverem(ii-1,Nb)+1;
	 cout_p = cout + M*(s_ii*Nblock-b2)+w*M*N ;
	 cbuf_p = plan.cbuf + M*(Nblock+b2)     +w*M*Nblocke;
	 for (int m=0; m<M; m++)
	 {
	    for (int n=0;n<b2;n++)
	    {
	       cout_p[m+n*M][0] += cbuf_p[m+n*M][0];
	       cout_p[m+n*M][1] += cbuf_p[m+n*M][1];
	    }
	 }
	 
      }
             
   }

   
}

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_done)(LTFAT_NAME(dgt_ola_plan) plan)
{
   LTFAT_NAME(dgt_long_done)(plan.plan);

   ltfat_free(plan.cbuf);
   ltfat_free(plan.gext);
   ltfat_free(plan.buf);

}



LTFAT_EXTERN LTFAT_NAME(dgtreal_ola_plan)
LTFAT_NAME(dgtreal_ola_init)(const LTFAT_REAL *g, const int gl,
			 const int W, const int a, const int M, const int bl,
			 unsigned flags)
{

   LTFAT_NAME(dgtreal_ola_plan) plan;

   plan.bl = bl;
   plan.gl = gl;
   plan.W  = W;
   const int M2=M/2+1;

   const int Lext    = bl+gl;
   const int Nblocke = Lext/a;
   
   plan.buf  = ltfat_malloc(Lext*W*sizeof(LTFAT_REAL));
   plan.gext = ltfat_malloc(Lext*sizeof(LTFAT_REAL));   
   plan.cbuf = ltfat_malloc(M2*Nblocke*W*sizeof(LTFAT_COMPLEX));
   
   LTFAT_NAME(fir2long_r)(g, gl, Lext, plan.gext);
   
   /* Zero the last part of the buffer, it will always be zero. */
   for (int w=0; w<W; w++)
   {      
      for (int jj=bl; jj<Lext;jj++)
      {
	 plan.buf[jj+w*Lext]=0.0;
      }
   }
   
   plan.plan =
      LTFAT_NAME(dgtreal_long_init)((const LTFAT_REAL*)plan.buf,
				(const LTFAT_REAL*)plan.gext,
				Lext, W, a, M,
				plan.cbuf, flags);
   
   return (plan);

}







LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_execute)(const LTFAT_NAME(dgtreal_ola_plan) plan,
			    const LTFAT_REAL *f, const int L,
			    LTFAT_COMPLEX *cout)

{
   const int bl      = plan.bl;
   const int gl      = plan.gl;
   const int a       = plan.plan.a;
   const int M       = plan.plan.M;
   const int N       = L/a;
   const int Lext    = bl+gl;
   const int Nb      = L/bl;
   const int b2      = gl/a/2;
   const int Nblock  = bl/a;
   const int Nblocke = Lext/a;
   const int W       = plan.W;
   const int M2      = M/2+1;

   /* Zero the output array, as we will be adding to it */
   for (int ii=0; ii<M2*N*W; ii++)
   {
      cout[ii][0]=0.0;
      cout[ii][1]=0.0;
   }

   
   for (int ii=0; ii<Nb; ii++)
   {
      int s_ii;

      /* Copy to working buffer. */
      for (int w=0; w<W; w++)
      {
	 memcpy(plan.buf+Lext*w,f+ii*bl+w*L,sizeof(LTFAT_REAL)*bl);
      }
      
      /* Execute the short DGTREAL */
      LTFAT_NAME(dgtreal_long_execute)(plan.plan);

      /* Place the results */
      for (int w=0; w<W; w++)
      {
	 LTFAT_COMPLEX *cout_p;
	 LTFAT_COMPLEX *cbuf_p;

	 /* Place large block */
	 cout_p = cout + ii*M2*Nblock+w*M2*N ;
	 cbuf_p = plan.cbuf +             w*M2*Nblocke;
	 for (int m=0; m<M2; m++)
	 {
	    for (int n=0;n<Nblock;n++)
	    {
	       cout_p[m+n*M2][0] += cbuf_p[m+n*M2][0];
	       cout_p[m+n*M2][1] += cbuf_p[m+n*M2][1];
	    }
	 }

	 /* Small block + */
	 s_ii=positiverem(ii+1,Nb);
	 cout_p = cout + s_ii*M2*Nblock+w*M2*N ;
	 cbuf_p = plan.cbuf +      M2*Nblock+w*M2*Nblocke;
	 for (int m=0; m<M2; m++)
	 {
	    for (int n=0;n<b2;n++)
	    {
	       cout_p[m+n*M2][0] += cbuf_p[m+n*M2][0];
	       cout_p[m+n*M2][1] += cbuf_p[m+n*M2][1];
	    }
	 }

	 
	 /* Small block - */
	 s_ii=positiverem(ii-1,Nb)+1;
	 cout_p = cout + M2*(s_ii*Nblock-b2)+w*M2*N ;
	 cbuf_p = plan.cbuf + M2*(     Nblock+b2)+w*M2*Nblocke;
	 for (int m=0; m<M2; m++)
	 {
	    for (int n=0;n<b2;n++)
	    {
	       cout_p[m+n*M2][0] += cbuf_p[m+n*M2][0];
	       cout_p[m+n*M2][1] += cbuf_p[m+n*M2][1];
	    }
	 }
	 
      }
             
   }
}


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_done)(LTFAT_NAME(dgtreal_ola_plan) plan)
{
   LTFAT_NAME(dgtreal_long_done)(plan.plan);

   ltfat_free(plan.cbuf);
   ltfat_free(plan.gext);
   ltfat_free(plan.buf);

}
