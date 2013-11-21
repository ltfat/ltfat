#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(ufilterbank_fft)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			    const int L, const int gl,
			    const int W, const int a, const int M,
			    LTFAT_COMPLEX *cout)
{

   /* ----- Initialization ------------ */

   const int N=L/a;

   LTFAT_COMPLEX *gwork = (LTFAT_COMPLEX*)ltfat_malloc(L*M*sizeof(LTFAT_COMPLEX));

   LTFAT_COMPLEX *work = (LTFAT_COMPLEX*)ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

   LTFAT_FFTW(plan) plan_g =
      LTFAT_FFTW(plan_many_dft)(1, &L, M,
				gwork, NULL,
				1, L,
				gwork, NULL,
				1, L,
				FFTW_FORWARD, FFTW_ESTIMATE);

      LTFAT_FFTW(plan_dft_1d)(L, gwork, gwork,
			      FFTW_FORWARD, FFTW_ESTIMATE);

   LTFAT_FFTW(plan) plan_w =
       LTFAT_FFTW(plan_dft_1d)(L, work, work,
			      FFTW_FORWARD, FFTW_ESTIMATE);

   LTFAT_FFTW(plan) plan_c =
      LTFAT_FFTW(plan_many_dft)(1, &N, M*W,
				cout, NULL,
				1, N,
				cout, NULL,
				1, N,
				FFTW_BACKWARD, FFTW_ESTIMATE);

   const LTFAT_REAL scalconst = 1.0/L;

   /* ----- Main -------------------------- */

   /* Extend g and copy to work buffer */
   for (int m=0; m<M; m++)
   {
      LTFAT_NAME(fir2long_c)(g+m*gl, gl, L, gwork+m*L);
   }

   LTFAT_FFTW(execute)(plan_g);

   for (int w=0; w<W; w++)
   {
      memcpy(work,f+L*w,sizeof(LTFAT_COMPLEX)*L);
      LTFAT_FFTW(execute)(plan_w);

      for (int m=0; m<M; m++)
      {
	 for (int n=0; n<N; n++)
	 {
        cout[n+m*N+w*N*M]=(LTFAT_COMPLEX) 0.0;
	    //cout[n+m*N+w*N*M][0]=0.0;
	    //cout[n+m*N+w*N*M][1]=0.0;
	    for (int k=0; k<a; k++)
	    {
	       const int l=n+k*N;
	       cout[n+m*N+w*N*M] += work[l]*gwork[l+m*L]*scalconst;
	       //const LTFAT_REAL tmp0 = work[l][0]*gwork[l+m*L][0]-work[l][1]*gwork[l+m*L][1];
	       //const LTFAT_REAL tmp1 = work[l][0]*gwork[l+m*L][1]+work[l][1]*gwork[l+m*L][0];
	       //cout[n+m*N+w*N*M][0]+=tmp0*scalconst;
	       //cout[n+m*N+w*N*M][1]+=tmp1*scalconst;
	    }
	 }
      }
   }


   LTFAT_FFTW(execute)(plan_c);



   ltfat_free(work);
   ltfat_free(gwork);

}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                        const size_t L, const size_t a,
                        LTFAT_COMPLEX *cout)
{
    const size_t Lc = L/a;
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(Lc, (LTFAT_COMPLEX*)cout, (LTFAT_COMPLEX*) cout,
                                                       FFTW_BACKWARD, FFTW_ESTIMATE);
    LTFAT_NAME(convsub_fft_plan)(F,G,L,a,cout, &plan_c);
    LTFAT_FFTW(destroy_plan)(plan_c);
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_plan)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                        const size_t L, const size_t a,
                        LTFAT_COMPLEX *cout, LTFAT_FFTW(plan)* p)
{
    const size_t Lc = L/a;
    const LTFAT_REAL scalconst = (LTFAT_REAL) 1.0/(L);
    LTFAT_COMPLEX* GPtrTmp = (LTFAT_COMPLEX*) G;
    LTFAT_COMPLEX* FPtrTmp = (LTFAT_COMPLEX*) F;

    memset(cout,0,Lc*sizeof(LTFAT_COMPLEX));

    for(size_t jj=0;jj<a;jj++)
    {
       for(size_t ii=0;ii<Lc;ii++)
       {
          cout[ii] += *GPtrTmp++**FPtrTmp++;
       }
    }

    for(size_t ii=0;ii<Lc;ii++)
    {
       cout[ii] *= scalconst;
    }

    LTFAT_FFTW(execute_dft)(*p,cout,cout);

}


LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl)(const LTFAT_COMPLEX *F, const size_t L,
                          const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                          const double a, const double realonly, LTFAT_COMPLEX *cout)
{


   const size_t Lc = (size_t) floor(L/a + 0.5);
   LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(Lc, cout, cout,
                                                       FFTW_BACKWARD, FFTW_ESTIMATE);

//LTFAT_FFTW(plan) plan_c = NULL;
   LTFAT_NAME(convsub_fftbl_plan)(F,L,G, Lg, foff, a, realonly, cout, &plan_c);

   LTFAT_FFTW(destroy_plan)(plan_c);

}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_plan)(const LTFAT_COMPLEX *F, const size_t L,
                          const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                          const double a, const double realonly, LTFAT_COMPLEX *cout, LTFAT_FFTW(plan)* p)
{

   const size_t Lc = (size_t) floor(L/a + 0.5);

   const size_t tmpLen = (size_t) ceil(Lg/((double)Lc))*Lc;
   const LTFAT_REAL scalconst = (LTFAT_REAL) 1.0/(L);
   LTFAT_COMPLEX *tmp = (LTFAT_COMPLEX*)ltfat_malloc(tmpLen*sizeof(LTFAT_COMPLEX));

   //memset(cout,0,Lc*sizeof(LTFAT_COMPLEX));
   memset(tmp,0,tmpLen*sizeof(LTFAT_COMPLEX));

   LTFAT_COMPLEX *tmpPtr = tmp;
   ptrdiff_t foffTmp = foff;
   size_t tmpLg = Lg;

   // Copy samples of F according to range of G
   if(foffTmp<0)
   {
       size_t toCopy = min_pt(-foffTmp,tmpLg);
       memcpy(tmpPtr,F+L+foffTmp,toCopy*sizeof(LTFAT_COMPLEX));
       tmpPtr+=toCopy;
       tmpLg-=toCopy;
       foffTmp = 0;
   }

   if(foffTmp+tmpLg>L)
   {
       ptrdiff_t over = foffTmp+tmpLg - L;
       memcpy(tmpPtr+Lg-over,F,over*sizeof(LTFAT_COMPLEX));
       tmpLg -=over;
   }

   memcpy(tmpPtr,F+foffTmp,tmpLg*sizeof(LTFAT_COMPLEX));

   // Do the filtering
   for(size_t ii=0;ii<Lg;ii++)
   {
      tmp[ii] *= G[ii];
   }

   // Do the folding
   for(size_t jj=1;jj<tmpLen/Lc;jj++)
   {
      for(size_t ii=0;ii<Lc;ii++)
      {
         tmp[ii] += tmp[jj*Lc+ii];
      }
   }

   LTFAT_NAME_COMPLEX(circshift)(tmp,cout,Lc,foff);


    for(size_t ii=0;ii<Lc;ii++)
    {
       cout[ii] *= scalconst;
    }

   // ifft
   LTFAT_FFTW(execute_dft)(*p,cout,cout);


   if(realonly>1e-3)
   {
      const ptrdiff_t foffconj = positiverem(L-foff-Lg,L)+1;
      LTFAT_COMPLEX *Gconj = (LTFAT_COMPLEX*)ltfat_malloc(Lg*sizeof(LTFAT_COMPLEX));
      for(size_t ii=0;ii<Lg;ii++)
      {
         Gconj[ii] = (LTFAT_COMPLEX) conj((double _Complex)G[Lg-1-ii]);
      }

      LTFAT_NAME(convsub_fftbl_plan)(F, L, Gconj, Lg, foffconj, a, 0, tmp, p);
      for(size_t ii=0;ii<Lc;ii++)
      {
         cout[ii] = (cout[ii] + tmp[ii])/2.0;
      }
      ltfat_free(Gconj);
   }

   ltfat_free(tmp);
}


// Inverse
LTFAT_EXTERN void
LTFAT_NAME(upconv_fft)(const LTFAT_COMPLEX *c, const size_t Lc,
                       const LTFAT_COMPLEX *G, const size_t a,
                       LTFAT_COMPLEX *Fout)
{
    LTFAT_COMPLEX* cbuffer = ltfat_malloc(Lc*sizeof(LTFAT_COMPLEX));
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(Lc, (LTFAT_COMPLEX*)cbuffer, (LTFAT_COMPLEX*) cbuffer,
                                                       FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_NAME(upconv_fft_plan)(c,Lc,G,a,Fout,&plan_c,cbuffer);
    LTFAT_FFTW(destroy_plan)(plan_c);
    ltfat_free(cbuffer);
}

LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_plan)(const LTFAT_COMPLEX *c, const size_t Lc,
                            const LTFAT_COMPLEX *G, const size_t a,
                            LTFAT_COMPLEX *Fout, LTFAT_FFTW(plan) *p,
                            LTFAT_COMPLEX *cbuffer
                            )
{
   memcpy(cbuffer,c,Lc*sizeof(LTFAT_COMPLEX));
   LTFAT_FFTW(execute_dft)(*p,cbuffer,cbuffer);

   LTFAT_COMPLEX* GPtrTmp = (LTFAT_COMPLEX*) G;
   LTFAT_COMPLEX* FPtrTmp = (LTFAT_COMPLEX*) Fout;

   for(size_t jj=0;jj<a;jj++)
   {
      for(size_t ii=0;ii<Lc;ii++)
      {
         // Really readable ;)
         *FPtrTmp++ += LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++)*cbuffer[ii];
      }
   }
}

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl)(const LTFAT_COMPLEX *c, const size_t Lc,
                            const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_COMPLEX *Fout)
{
    LTFAT_COMPLEX* cbuffer = ltfat_malloc(Lc*sizeof(LTFAT_COMPLEX));
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(Lc, (LTFAT_COMPLEX*)cbuffer, (LTFAT_COMPLEX*) cbuffer,
                                                       FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_NAME(upconv_fftbl_plan)(c,Lc,G,Lg,foff,a,realonly,Fout,&plan_c,cbuffer);
    LTFAT_FFTW(destroy_plan)(plan_c);
    ltfat_free(cbuffer);
}

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_plan)(const LTFAT_COMPLEX *c, const size_t Lc,
                            const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_COMPLEX *Fout,
                            LTFAT_FFTW(plan) *p,LTFAT_COMPLEX *cbuffer)
{
   memcpy(cbuffer,c,Lc*sizeof(LTFAT_COMPLEX));
   LTFAT_FFTW(execute_dft)(*p,cbuffer,cbuffer);

   size_t L = (size_t) floor(Lc*a + 0.5); // round

   LTFAT_NAME_COMPLEX(circshift)(cbuffer,cbuffer,Lc,-foff);

   LTFAT_COMPLEX* GPtrTmp = (LTFAT_COMPLEX*) G;
   LTFAT_COMPLEX* FPtrTmp = (LTFAT_COMPLEX*) Fout;
   LTFAT_COMPLEX* CPtrTmp = (LTFAT_COMPLEX*) cbuffer;

   // Determine range of G
   ptrdiff_t foffTmp = foff;
   ptrdiff_t tmpLg = Lc<Lg?Lc:Lg;
   ptrdiff_t over = 0;
   if(foffTmp+tmpLg>(ptrdiff_t)L)
   {
      over = foffTmp+tmpLg - (ptrdiff_t)L;
   }


   if(foffTmp<0)
   {
      size_t toCopy = (-foffTmp)<tmpLg?-foffTmp:tmpLg;
      FPtrTmp = Fout+L+foffTmp;
      for(size_t ii=0;ii<toCopy;ii++)
      {
         FPtrTmp[ii]+=*CPtrTmp++ * LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++);
      }

      tmpLg-=toCopy;
      FPtrTmp = Fout;
      foffTmp = 0;
   }

   FPtrTmp = (LTFAT_COMPLEX*) Fout+foffTmp;
   for(ptrdiff_t ii=0;ii<tmpLg-over;ii++)
   {
      FPtrTmp[ii]+=*CPtrTmp++ * LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++);
   }

   for(ptrdiff_t ii=0;ii<over;ii++)
   {
      Fout[ii]+=*CPtrTmp++ * LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++);
   }


   if(realonly>1e-3)
   {
      const ptrdiff_t foffconj = positiverem(L-foff-Lg,L)+1;
      LTFAT_COMPLEX *Gconj = (LTFAT_COMPLEX*)ltfat_malloc(Lg*sizeof(LTFAT_COMPLEX));
      LTFAT_NAME_COMPLEX(reverse_array)((LTFAT_COMPLEX *)G,Gconj,Lg);
      LTFAT_NAME_COMPLEX(conjugate_array)(Gconj,Gconj,Lg);

      LTFAT_NAME(upconv_fftbl_plan)(c, Lc, Gconj, Lg, foffconj, a, 0, Fout, p, cbuffer);
      ltfat_free(Gconj);
   }

}
