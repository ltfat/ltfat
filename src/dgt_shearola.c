#include "config.h"
#include <string.h>
#include <complex.h>
#include "fftw3.h"
#include "ltfat.h"


LTFAT_EXTERN LTFAT_NAME(dgt_shearola_plan)
LTFAT_NAME(dgt_shearola_init)(const LTFAT_COMPLEX *g, const int gl,
			      const int W, const int a, const int M, 
			      const int s0, const int s1, const int br,
			      const int bl,
			      unsigned flags)
{

   LTFAT_NAME(dgt_shearola_plan) plan;

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
      LTFAT_NAME(dgt_shear_init)((const LTFAT_COMPLEX*)plan.buf,
				 (const LTFAT_COMPLEX*)plan.gext,
				 Lext, W, a, M,
				 s0, s1, br,
				 plan.cbuf, flags);
   
   return (plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_execute)(const LTFAT_NAME(dgt_shearola_plan) plan,
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
      LTFAT_NAME(dgt_shear_execute)(plan.plan);

      /* Place the results */
      for (int w=0; w<W; w++)
      {
	 /* Place large block */
	 LTFAT_COMPLEX *cout_p = cout +      ii*M*Nblock+w*M*N ;
	 LTFAT_COMPLEX *cbuf_p = plan.cbuf +  w*M*Nblocke;
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
LTFAT_NAME(dgt_shearola_done)(LTFAT_NAME(dgt_shearola_plan) plan)
{
   LTFAT_NAME(dgt_shear_done)(plan.plan);

   /* ltfat_free(plan.cbuf); */
   ltfat_free(plan.gext);
   ltfat_free(plan.buf);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			 const int L, const int gl, const int W, const int a, const int M,
			 const int s0, const int s1, const int br, const int bl,
			 LTFAT_COMPLEX *cout)
{ 

   LTFAT_NAME(dgt_shearola_plan) plan = LTFAT_NAME(dgt_shearola_init)(
      g,gl,W,a,M,s0,s1,br,bl,FFTW_ESTIMATE);

   LTFAT_NAME(dgt_shearola_execute)(plan,f,L,cout);

   LTFAT_NAME(dgt_shearola_done)(plan);

}
