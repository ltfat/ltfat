#include "ltfat.h"
#include "ltfat_types.h"


/**
* FFT filterbank routines
*/

struct LTFAT_NAME(convsub_fft_plan_struct)
{
    const ltfatInt L;
    const ltfatInt W;
    const ltfatInt a;
    const LTFAT_FFTW(plan) p_c;
};

struct LTFAT_NAME(convsub_fftbl_plan_struct)
{
    const ltfatInt L;
    const ltfatInt Gl;
    const ltfatInt W;
    const double a;
    const LTFAT_FFTW(plan) p_c;
    LTFAT_COMPLEX* buf;
    const ltfatInt bufLen;
};

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt W, const ltfatInt a[], const ltfatInt M,
                           LTFAT_COMPLEX *cout[])
{
    for(ltfatInt m =0; m<M; m++)
    {
        LTFAT_NAME(convsub_fft)(F,G[m],L,W,a[m],cout[m]);
    }
}


LTFAT_EXTERN void
LTFAT_NAME(filterbank_fft_execute)(LTFAT_NAME(convsub_fft_plan) p[],
                                   const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                                   const ltfatInt M, LTFAT_COMPLEX *cout[])
{

    for(ltfatInt m =0; m<M; m++)
    {
        LTFAT_NAME(convsub_fft_execute)(p[m],F,G[m],cout[m]);
    }
}


LTFAT_EXTERN LTFAT_NAME(convsub_fft_plan)
LTFAT_NAME(convsub_fft_init)(const ltfatInt L, const ltfatInt W,
                             const ltfatInt a, const LTFAT_COMPLEX *cout)
{
    const ltfatInt N = L/a;

    int Nint = (int) N;

    const LTFAT_FFTW(iodim) dims = {.n = Nint, .is = 1, .os = 1};
    const LTFAT_FFTW(iodim) howmany_dims = {.n = W,.is = Nint, .os = Nint};

    LTFAT_COMPLEX* coutNc = (LTFAT_COMPLEX*) cout;
    LTFAT_FFTW(plan) p_many =
        LTFAT_FFTW(plan_guru_dft)(1, &dims, 1, &howmany_dims,
                                  coutNc, coutNc,
                                  FFTW_BACKWARD, FFTW_ESTIMATE);


    struct LTFAT_NAME(convsub_fft_plan_struct) p_struct =
    { .L = L, .a = a, .W = W, .p_c = p_many };

    LTFAT_NAME(convsub_fft_plan) p = ltfat_malloc(sizeof*p);
    memcpy(p,&p_struct,sizeof*p);
    return p;
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_done)(LTFAT_NAME(convsub_fft_plan) p)
{
    LTFAT_FFTW(destroy_plan)(p->p_c);
    ltfat_free(p);
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_execute)(const LTFAT_NAME(convsub_fft_plan) p,
                                const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                                LTFAT_COMPLEX *cout)
{
    const ltfatInt L = p->L;
    const ltfatInt W = p->W;
    const ltfatInt a = p->a;
    const ltfatInt N = L/a;
    const LTFAT_REAL scalconst = (LTFAT_REAL) (1.0/L);

    memset(cout,0,W*N*sizeof*cout);

    for(ltfatInt w=0; w<W; w++)
    {
        LTFAT_COMPLEX* GPtrTmp = (LTFAT_COMPLEX*) G;
        LTFAT_COMPLEX* FPtrTmp = (LTFAT_COMPLEX*) (F + w*L);
        for(ltfatInt jj=0; jj<a; jj++)
        {
            for(ltfatInt n=0; n<N; n++)
            {
                cout[w*N+n] += *GPtrTmp++**FPtrTmp++;
            }
        }
    }

    for(ltfatInt ii=0; ii<N*W; ii++)
    {
        cout[ii] *= scalconst;
    }

    LTFAT_FFTW(execute_dft)(p->p_c,cout,cout);

}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                        const ltfatInt L, const ltfatInt W,
                        const ltfatInt a, LTFAT_COMPLEX *cout)
{
    LTFAT_NAME(convsub_fft_plan) p = LTFAT_NAME(convsub_fft_init)(L,W,a,cout);
    LTFAT_NAME(convsub_fft_execute)(p,F,G,cout);
    LTFAT_NAME(convsub_fft_done)(p);
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
        LTFAT_NAME(convsub_fftbl)(F,G[m],L,Gl[m],W,a[m],
                                  foff[m],realonly[m],cout[m]);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fftbl_execute)(LTFAT_NAME(convsub_fftbl_plan) p[],
                                     const LTFAT_COMPLEX *F,
                                     const LTFAT_COMPLEX *G[],
                                     const ltfatInt M, const ltfatInt foff[],
                                     const int realonly[], LTFAT_COMPLEX *cout[])
{
    for(ltfatInt m =0; m<M; m++)
    {
        LTFAT_NAME(convsub_fftbl_execute)(p[m],F,G[m],foff[m],realonly[m],cout[m]);
    }

}


LTFAT_EXTERN LTFAT_NAME(convsub_fftbl_plan)
LTFAT_NAME(convsub_fftbl_init)( const ltfatInt L, const ltfatInt Gl,
                                const ltfatInt W, const double a,
                                const LTFAT_COMPLEX *cout)
{
    const ltfatInt N = (ltfatInt) floor(L/a + 0.5);

    int Nint = (int) N;

    const LTFAT_FFTW(iodim) dims = {.n = Nint, .is = 1, .os = 1};
    const LTFAT_FFTW(iodim) howmany_dims = {.n = W,.is = Nint, .os = Nint};

   LTFAT_COMPLEX *coutNc = (LTFAT_COMPLEX *) cout;
    LTFAT_FFTW(plan) p_many =
        LTFAT_FFTW(plan_guru_dft)(1, &dims, 1, &howmany_dims,
                                  coutNc, coutNc,
                                  FFTW_BACKWARD, FFTW_ESTIMATE);

    const ltfatInt bufLen = (ltfatInt) ceil(Gl/((double)N))*N;

    LTFAT_COMPLEX* buf = NULL;
    if(bufLen)
       buf = ltfat_malloc(bufLen*sizeof*buf);

    struct LTFAT_NAME(convsub_fftbl_plan_struct) p_struct =
    {
        .L = L, .Gl = Gl, .a = a, .W = W, .p_c = p_many , .bufLen = bufLen,
        .buf = buf
    };

    LTFAT_NAME(convsub_fftbl_plan) p = ltfat_malloc(sizeof*p);
    memcpy(p,&p_struct,sizeof*p);
    return p;
}

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_done)( LTFAT_NAME(convsub_fftbl_plan) p)
{
    LTFAT_FFTW(destroy_plan)(p->p_c);
    if(p->buf!=NULL) ltfat_free(p->buf);
    ltfat_free(p);
}


LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_execute)(const LTFAT_NAME(convsub_fftbl_plan) p,
                                  const LTFAT_COMPLEX *F,
                                  const LTFAT_COMPLEX *G,
                                  const ltfatInt foff,
                                  const int realonly,
                                  LTFAT_COMPLEX *cout)
{

    const ltfatInt L = p->L;
    const ltfatInt Gl = p->Gl;
    const ltfatInt W = p->W;
    const double a = p->a;
    // Output length
    const ltfatInt N = (ltfatInt) floor(L/a + 0.5);
    // Bail out in degenerate case
    if(!Gl)
    {
        memset(cout,0,W*N*sizeof*cout);
        return;
    }

    LTFAT_COMPLEX* tmp = p->buf;
    const ltfatInt tmpLen = p->bufLen;
    //const ltfatInt tmpLen = (ltfatInt) ceil(Gl/((double)N))*N;
    const LTFAT_REAL scalconst = (LTFAT_REAL) 1.0/(L);
    // LTFAT_COMPLEX *tmp = ltfat_calloc(tmpLen,sizeof*tmp);


    for(ltfatInt w=0; w<W; w++)
    {
        // First Gl elements of tmp is copied from F so,
        // zero only the part which wont be written to.
        memset(tmp+Gl,0,(tmpLen-Gl)*sizeof*tmp);
        LTFAT_COMPLEX *tmpPtr = tmp;
        ltfatInt foffTmp = foff;
        ltfatInt tmpLg = Gl;

        // Copy samples of F according to range of G
        if(foffTmp<0)
        {
            ltfatInt toCopy = imin(-foffTmp,tmpLg);
            memcpy(tmpPtr,F+(w+1)*L+foffTmp,toCopy*sizeof*F);
            tmpPtr+=toCopy;
            tmpLg-=toCopy;
            foffTmp = 0;
        }

        if(foffTmp+tmpLg>L)
        {
            ltfatInt over = foffTmp+tmpLg - L;
            memcpy(tmpPtr+Gl-over,F+w*L,over*sizeof*F);
            tmpLg -=over;
        }

        memcpy(tmpPtr,F+w*L+foffTmp,tmpLg*sizeof*F);

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

        // Do the circshift
         LTFAT_NAME_COMPLEX(circshift)(tmp,cout+w*N,N,foff);
        //LTFAT_NAME_COMPLEX(circshift)(tmp,cout+w*N,N,-Gl/2);
        // memcpy(cout+w*N,tmp,N*sizeof*cout);
    }

    for(ltfatInt ii=0; ii<W*N; ii++)
    {
        cout[ii] *= scalconst;
    }

    // ifft
    LTFAT_FFTW(execute_dft)(p->p_c,cout,cout);


    if(realonly)
    {
        // Involute the filter and call the function again
        const ltfatInt foffconj = -L + positiverem(L-foff-Gl,L)+1;
        LTFAT_COMPLEX *Gconj = ltfat_malloc(Gl*sizeof*Gconj);
        LTFAT_COMPLEX *cout2 = ltfat_malloc(W*N*sizeof*cout2);
        for(ltfatInt ii=0; ii<Gl; ii++)
        {
            Gconj[ii] = (LTFAT_COMPLEX) conj((double _Complex)G[Gl-1-ii]);
        }

        LTFAT_NAME(convsub_fftbl_execute)(p, F, Gconj, foffconj, 0, cout2);

        // Scale
        for(ltfatInt ii=0; ii<W*N; ii++)
        {
            cout[ii] = (cout[ii] + cout2[ii])/2.0;
        }
        ltfat_free(Gconj);
        ltfat_free(cout2);
    }
}


LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl)(const LTFAT_COMPLEX *F,  const LTFAT_COMPLEX *G,
                          const ltfatInt L, const ltfatInt Gl, const ltfatInt W,
                          const double a, const ltfatInt foff,
                          const int realonly, LTFAT_COMPLEX *cout)
{
    LTFAT_NAME(convsub_fftbl_plan) p =
        LTFAT_NAME(convsub_fftbl_init)( L, Gl, W, a, cout);

    LTFAT_NAME(convsub_fftbl_execute)(p, F, G, foff, realonly, cout);

    LTFAT_NAME(convsub_fftbl_done)(p);
}


LTFAT_EXTERN void
LTFAT_NAME(ufilterbank_fft)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                            const ltfatInt L, const ltfatInt Gl,
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
        LTFAT_NAME(fir2long_c)(g+m*Gl, Gl, L, gwork+m*L);
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
