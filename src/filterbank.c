#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"


/**
* FFT filterbank routines
*/

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt W, const ltfatInt a[], const ltfatInt M,
                           LTFAT_COMPLEX *cout[])
{
    for(ltfatInt m =0; m<M; m++)
    {
        ltfatInt N = L/a[m];
        // First col of cPtrs[m] is used as a buffer to assure correct memory aligment
        for(ltfatInt w =1; w<W; w++)
        {
            LTFAT_NAME(convsub_fft)(F+w*L,G[m],L,a[m],cout[m]);
            memcpy(cout[m] + w*N, cout[m], N*sizeof(LTFAT_COMPLEX));
        }

        LTFAT_NAME(convsub_fft)(F,G[m],L,a[m],cout[m]);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fft_plans)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                                 const ltfatInt L, const ltfatInt W, const ltfatInt a[], const ltfatInt M,
                                 LTFAT_COMPLEX *cout[], LTFAT_FFTW(plan) p[])
{
    for(ltfatInt m =0; m<M; m++)
    {
        ltfatInt N = L/a[m];
        // First col of cPtrs[m] is used as a buffer to assure correct memory aligment
        for(ltfatInt w =1; w<W; w++)
        {
            LTFAT_NAME(convsub_fft_plan)(F+w*L,G[m],L,a[m],cout[m],&p[m]);
            memcpy(cout[m] + w*N, cout[m], N*sizeof(LTFAT_COMPLEX));
        }

        LTFAT_NAME(convsub_fft_plan)(F,G[m],L,a[m],cout[m],&p[m]);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                        const ltfatInt L, const ltfatInt a, LTFAT_COMPLEX *cout)
{
    const ltfatInt N = L/a;
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(N, (LTFAT_COMPLEX*)cout, (LTFAT_COMPLEX*) cout,
                               FFTW_BACKWARD, FFTW_ESTIMATE);
    LTFAT_NAME(convsub_fft_plan)(F,G,L,a,cout, &plan_c);
    LTFAT_FFTW(destroy_plan)(plan_c);
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_plan)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                             const ltfatInt L, const ltfatInt a, LTFAT_COMPLEX *cout,
                             LTFAT_FFTW(plan) *p)
{
    const ltfatInt N = L/a;
    const LTFAT_REAL scalconst = (LTFAT_REAL) 1.0/(L);
    LTFAT_COMPLEX* GPtrTmp = (LTFAT_COMPLEX*) G;
    LTFAT_COMPLEX* FPtrTmp = (LTFAT_COMPLEX*) F;

    memset(cout,0,N*sizeof(LTFAT_COMPLEX));

    for(ltfatInt jj=0; jj<a; jj++)
    {
        for(ltfatInt ii=0; ii<N; ii++)
        {
            cout[ii] += *GPtrTmp++**FPtrTmp++;
        }
    }

    for(ltfatInt ii=0; ii<N; ii++)
    {
        cout[ii] *= scalconst;
    }

    LTFAT_FFTW(execute_dft)(*p,cout,cout);

}

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fftbl)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                             const ltfatInt L, const ltfatInt Gl[],
                             const ltfatInt W, const double a[], const ltfatInt M,
                             const ltfatInt foff[], const int realonly[],
                             LTFAT_COMPLEX *cout[])
{


    for(ltfatInt m =0; m<M; m++)
    {
        const ltfatInt N = (ltfatInt) floor(L/a[m] + 0.5);
        for(ltfatInt w =1; w<W; w++)
        {
            // Using the first col of c as a temp array.
            LTFAT_NAME(convsub_fftbl)(F+w*L,G[m],L,Gl[m],a[m],foff[m],realonly[m],cout[m]);

            // Copy to an appropriate position
            memcpy(cout[m] + w*N,cout[m],N*sizeof*cout[m]);
        }

        // Working with the first col only.
        LTFAT_NAME(convsub_fftbl)(F,G[m],L,Gl[m],a[m],foff[m],realonly[m],cout[m]);
    }

}

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fftbl_plans)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                                   const ltfatInt L, const ltfatInt Gl[],
                                   const ltfatInt W, const double a[], const ltfatInt M,
                                   const ltfatInt foff[], const int realonly[],
                                   LTFAT_COMPLEX *cout[],LTFAT_FFTW(plan) p[])
{

    for(ltfatInt m =0; m<M; m++)
    {
        const ltfatInt N = (ltfatInt) floor(L/a[m] + 0.5);
        for(ltfatInt w =1; w<W; w++)
        {
            // Using the first col of c as a temp array.
            LTFAT_NAME(convsub_fftbl_plan)(F+w*L,G[m],L,Gl[m],a[m],foff[m],realonly[m],cout[m],&p[m]);

            // Copy to an appropriate position
            memcpy(cout[m] + w*N,cout[m],N*sizeof*cout[m]);
        }

        // Working with the first col only.
        LTFAT_NAME(convsub_fftbl_plan)(F,G[m],L,Gl[m],a[m],foff[m],realonly[m],cout[m],&p[m]);
    }

}



LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl)(const LTFAT_COMPLEX *F,  const LTFAT_COMPLEX *G,
                          const ltfatInt L, const ltfatInt Gl, const double a,
                          const ltfatInt foff, const int realonly, LTFAT_COMPLEX *cout)
{
    const ltfatInt N = (ltfatInt) floor(L/a + 0.5);
    LTFAT_FFTW(plan) plan_c =  LTFAT_FFTW(plan_dft_1d)(N, cout, cout,
                               FFTW_BACKWARD, FFTW_ESTIMATE);

    LTFAT_NAME(convsub_fftbl_plan)(F,G,L, Gl, a, foff, realonly, cout, &plan_c);

    LTFAT_FFTW(destroy_plan)(plan_c);
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_plan)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G, const ltfatInt L,
                               const ltfatInt Gl, const double a, const ltfatInt foff,
                               const int realonly, LTFAT_COMPLEX *cout, LTFAT_FFTW(plan)* p)
{
    // Output length
    const ltfatInt N = (ltfatInt) floor(L/a + 0.5);

    const ltfatInt tmpLen = (ltfatInt) ceil(Gl/((double)N))*N;
    const LTFAT_REAL scalconst = (LTFAT_REAL) 1.0/(L);
    LTFAT_COMPLEX *tmp = ltfat_calloc(tmpLen,sizeof*tmp);

    LTFAT_COMPLEX *tmpPtr = tmp;
    ltfatInt foffTmp = foff;
    ltfatInt tmpLg = Gl;

    // Copy samples of F according to range of G
    if(foffTmp<0)
    {
        ltfatInt toCopy = imin(-foffTmp,tmpLg);
        memcpy(tmpPtr,F+L+foffTmp,toCopy*sizeof(LTFAT_COMPLEX));
        tmpPtr+=toCopy;
        tmpLg-=toCopy;
        foffTmp = 0;
    }

    if(foffTmp+tmpLg>L)
    {
        ltfatInt over = foffTmp+tmpLg - L;
        memcpy(tmpPtr+Gl-over,F,over*sizeof(LTFAT_COMPLEX));
        tmpLg -=over;
    }

    memcpy(tmpPtr,F+foffTmp,tmpLg*sizeof(LTFAT_COMPLEX));

    // Do the filtering
    for(ltfatInt ii=0; ii<Gl; ii++)
    {
        tmp[ii] *= G[ii];
    }

    // Do the folding
    for(ltfatInt jj=1; jj<tmpLen/N; jj++)
    {
        for(ltfatInt ii=0; ii<N; ii++)
        {
            tmp[ii] += tmp[jj*N+ii];
        }
    }

    LTFAT_NAME_COMPLEX(circshift)(tmp,cout,N,foff);


    for(ltfatInt ii=0; ii<N; ii++)
    {
        cout[ii] *= scalconst;
    }

    // ifft
    LTFAT_FFTW(execute_dft)(*p,cout,cout);


    if(realonly)
    {
        // Involute the filter and call the function again
        const ltfatInt foffconj = positiverem(L-foff-Gl,L)+1;
        LTFAT_COMPLEX *Gconj = ltfat_malloc(Gl*sizeof*Gconj);
        for(ltfatInt ii=0; ii<Gl; ii++)
        {
            Gconj[ii] = (LTFAT_COMPLEX) conj((double _Complex)G[Gl-1-ii]);
        }

        LTFAT_NAME(convsub_fftbl_plan)(F, Gconj, L, Gl, a, foffconj, false, tmp, p);

        // Scale
        for(ltfatInt ii=0; ii<N; ii++)
        {
            cout[ii] = (cout[ii] + tmp[ii])/2.0;
        }
        ltfat_free(Gconj);
    }

    ltfat_free(tmp);
}

LTFAT_EXTERN void
LTFAT_NAME(ufilterbank_fft)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                            const ltfatInt L, const ltfatInt gl,
                            const ltfatInt W, const ltfatInt a, const ltfatInt M,
                            LTFAT_COMPLEX *cout)
{

    /* ----- Initialization ------------ */

    const ltfatInt N=L/a;

    /* Downcasting to ints */
    int Lint = (int) L;
    int Nint = (int) N;

    LTFAT_COMPLEX *gwork = (LTFAT_COMPLEX*)ltfat_malloc(L*M*sizeof(LTFAT_COMPLEX));

    LTFAT_COMPLEX *work = (LTFAT_COMPLEX*)ltfat_malloc(L*sizeof(LTFAT_COMPLEX));

    LTFAT_FFTW(plan) plan_g =
        LTFAT_FFTW(plan_many_dft)(1, &Lint, M,
                                  gwork, NULL,
                                  1, Lint,
                                  gwork, NULL,
                                  1, Lint,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_FFTW(plan_dft_1d)(L, gwork, gwork,
                            FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_FFTW(plan) plan_w =
        LTFAT_FFTW(plan_dft_1d)(L, work, work,
                                FFTW_FORWARD, FFTW_ESTIMATE);

    LTFAT_FFTW(plan) plan_c =
        LTFAT_FFTW(plan_many_dft)(1, &Nint, M*W,
                                  cout, NULL,
                                  1, Nint,
                                  cout, NULL,
                                  1, Nint,
                                  FFTW_BACKWARD, FFTW_ESTIMATE);

    const LTFAT_REAL scalconst = 1.0/L;

    /* ----- Main -------------------------- */

    /* Extend g and copy to work buffer */
    for (ltfatInt m=0; m<M; m++)
    {
        LTFAT_NAME(fir2long_c)(g+m*gl, gl, L, gwork+m*L);
    }

    LTFAT_FFTW(execute)(plan_g);

    for (ltfatInt w=0; w<W; w++)
    {
        memcpy(work,f+L*w,sizeof(LTFAT_COMPLEX)*L);
        LTFAT_FFTW(execute)(plan_w);

        for (ltfatInt m=0; m<M; m++)
        {
            for (ltfatInt n=0; n<N; n++)
            {
                cout[n+m*N+w*N*M]=(LTFAT_COMPLEX) 0.0;

                for (ltfatInt k=0; k<a; k++)
                {
                    const ltfatInt l=n+k*N;
                    cout[n+m*N+w*N*M] += work[l]*gwork[l+m*L]*scalconst;
                }
            }
        }
    }


    LTFAT_FFTW(execute)(plan_c);


    LTFAT_SAFEFREEALL(work,gwork);
}
