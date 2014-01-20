#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt W, const ltfatInt a[], const ltfatInt M,
                           LTFAT_COMPLEX *F)
{
   // This is necessary since F us used as an accumulator
   memset(F,0,L*sizeof*F);

   for(ltfatInt m =0; m<M; m++)
    {
        ltfatInt N = L/a[m];
        for(ltfatInt w =0; w<W; w++)
        {
            LTFAT_NAME(upconv_fft)(cin[m]+w*N,G[m],L,a[m],F+w*L);
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft_plans)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt W, const ltfatInt a[], const ltfatInt M,
                           LTFAT_COMPLEX *F, LTFAT_FFTW(plan) p[], LTFAT_COMPLEX *cbuf[])
{
   // This is necessary since F us used as an accumulator
   memset(F,0,L*sizeof*F);

   for(ltfatInt m =0; m<M; m++)
    {
        ltfatInt N = L/a[m];
        for(ltfatInt w =0; w<W; w++)
        {
            LTFAT_NAME(upconv_fft_plan)(cin[m]+w*N,G[m],L,a[m],F+w*L,&p[m],cbuf[m]);
        }
    }
}


// Inverse
LTFAT_EXTERN void
LTFAT_NAME(upconv_fft)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                       const ltfatInt L, const ltfatInt a, LTFAT_COMPLEX *F)
{
    ltfatInt N = L/a;
    LTFAT_COMPLEX* cbuf = ltfat_malloc(N*sizeof*cbuf);
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(N, cbuf, cbuf,
                               FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_NAME(upconv_fft_plan)(cin,G,L,a,F,&plan_c,cbuf);

    LTFAT_FFTW(destroy_plan)(plan_c);
    ltfat_free(cbuf);
}

LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_plan)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                            const ltfatInt L, const ltfatInt a,
                            LTFAT_COMPLEX *F, LTFAT_FFTW(plan) *p,
                            LTFAT_COMPLEX *cbuf
                           )
{
    ltfatInt N = L/a;
    memcpy(cbuf,cin,N*sizeof*cin);
    // New array execution, inplace
    LTFAT_FFTW(execute_dft)(*p,cbuf,cbuf);

    for(ltfatInt jj=0; jj<a; jj++)
    {
        for(ltfatInt ii=0; ii<N; ii++)
        {
            // Really readable ;)
            *F++ += LTFAT_COMPLEXH_NAME(conj)(*G++)*cbuf[ii];
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fftbl)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt Gl[], const ltfatInt W, const double a[], const ltfatInt M,
                           const ltfatInt foff[], const int realonly[],
                           LTFAT_COMPLEX *F)
{
   // This is necessary since F us used as an accumulator
   memset(F,0,L*sizeof*F);

   for(ltfatInt m =0; m<M; m++)
    {
        ltfatInt N = (ltfatInt) floor(L/a[m] + 0.5);
        for(ltfatInt w =0; w<W; w++)
        {
            LTFAT_NAME(upconv_fftbl)(cin[m]+w*N,G[m],L,Gl[m],a[m],foff[m],realonly[m],F+w*L);
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fftbl_plans)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt Gl[], const ltfatInt W, const double a[], const ltfatInt M,
                           const ltfatInt foff[], const int realonly[],
                           LTFAT_COMPLEX *F, LTFAT_FFTW(plan) p[], LTFAT_COMPLEX *cbuf[])
{
   // This is necessary since F us used as an accumulator
   memset(F,0,L*sizeof*F);

   for(ltfatInt m =0; m<M; m++)
    {
        ltfatInt N = (ltfatInt) floor(L/a[m] + 0.5);
        for(ltfatInt w =0; w<W; w++)
        {
            LTFAT_NAME(upconv_fftbl_plan)(cin[m]+w*N,G[m],L,Gl[m],a[m],foff[m],realonly[m],F+w*L,&p[m],cbuf[m]);
        }
    }
}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                         const ltfatInt L, const ltfatInt Gl, const double a,
                         const ltfatInt foff, const int realonly,
                         LTFAT_COMPLEX *F)
{
    ltfatInt N = (ltfatInt) floor(L/a + 0.5);
    LTFAT_COMPLEX* cbuf = ltfat_malloc(N*sizeof*cbuf);
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(N, cbuf, cbuf,
                               FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_NAME(upconv_fftbl_plan)(cin,G,L,Gl,a,foff,realonly,F,&plan_c,cbuf);
    LTFAT_FFTW(destroy_plan)(plan_c);
    ltfat_free(cbuf);
}

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_plan)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                              const ltfatInt L, const ltfatInt Gl, const double a,
                              const ltfatInt foff, const int realonly, LTFAT_COMPLEX *F,
                              LTFAT_FFTW(plan) *p, LTFAT_COMPLEX *cbuf)
{
    ltfatInt N = (ltfatInt) floor(L/a + 0.5);
    memcpy(cbuf,cin,N*sizeof*cin);
    LTFAT_FFTW(execute_dft)(*p,cbuf,cbuf);

    LTFAT_NAME_COMPLEX(circshift)(cbuf,cbuf,N,-foff);

    const LTFAT_COMPLEX* GPtrTmp = G;
    LTFAT_COMPLEX* FPtrTmp = F;
    LTFAT_COMPLEX* CPtrTmp = cbuf;

    // Determine range of G
    ltfatInt foffTmp = foff;
    ltfatInt tmpLg = N<Gl?N:Gl;
    ltfatInt over = 0;
    if(foffTmp+tmpLg>(ltfatInt)L)
    {
        over = foffTmp+tmpLg - (ltfatInt)L;
    }


    if(foffTmp<0)
    {
        ltfatInt toCopy = (-foffTmp)<tmpLg?-foffTmp:tmpLg;
        FPtrTmp = F+L+foffTmp;
        for(ltfatInt ii=0; ii<toCopy; ii++)
        {
            FPtrTmp[ii]+=*CPtrTmp++ * LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++);
        }

        tmpLg-=toCopy;
        foffTmp = 0;
    }

    FPtrTmp = F+foffTmp;
    for(ltfatInt ii=0; ii<tmpLg-over; ii++)
    {
        FPtrTmp[ii]+=*CPtrTmp++ * LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++);
    }

    for(ltfatInt ii=0; ii<over; ii++)
    {
        F[ii]+=*CPtrTmp++ * LTFAT_COMPLEXH_NAME(conj)(*GPtrTmp++);
    }


    if(realonly)
    {
        const ltfatInt foffconj = positiverem(L-foff-Gl,L)+1;
        LTFAT_COMPLEX *Gconj = ltfat_malloc(Gl*sizeof*Gconj);
        LTFAT_NAME_COMPLEX(reverse_array)((LTFAT_COMPLEX *)G,Gconj,Gl);
        LTFAT_NAME_COMPLEX(conjugate_array)(Gconj,Gconj,Gl);

        LTFAT_NAME(upconv_fftbl_plan)(cin, Gconj, L, Gl, a, foffconj, false, F, p, cbuf);
        ltfat_free(Gconj);
    }

}

