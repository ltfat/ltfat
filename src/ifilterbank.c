#include "ltfat.h"
#include "ltfat_types.h"

struct LTFAT_NAME(upconv_fft_plan_struct)
{
    const ltfatInt L;
    const ltfatInt W;
    const ltfatInt a;
    const LTFAT_FFTW(plan) p_c;
    LTFAT_COMPLEX* buf;
    const ltfatInt bufLen;
};

struct LTFAT_NAME(upconv_fftbl_plan_struct)
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
LTFAT_NAME(ifilterbank_fft)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                            const ltfatInt L, const ltfatInt W, const ltfatInt a[],
                            const ltfatInt M, LTFAT_COMPLEX *F)
{
    // This is necessary since F us used as an accumulator
    memset(F, 0, L * W * sizeof * F);

    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_NAME(upconv_fft)(cin[m], G[m], L, W, a[m], F);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft_execute)(LTFAT_NAME(upconv_fft_plan) p[],
                                    const LTFAT_COMPLEX *cin[],
                                    const LTFAT_COMPLEX *G[],
                                    const ltfatInt M,
                                    LTFAT_COMPLEX *F )
{
    ltfatInt L = p[0]->L;
    ltfatInt W = p[0]->W;
    // This is necessary since F us used as an accumulator
    memset(F, 0, W * L * sizeof * F);

    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_NAME(upconv_fft_execute)(p[m], cin[m], G[m], F);
    }
}


// Inverse
LTFAT_EXTERN void
LTFAT_NAME(upconv_fft)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                       const ltfatInt L, const ltfatInt W, const ltfatInt a,
                       LTFAT_COMPLEX *F)
{
    LTFAT_NAME(upconv_fft_plan) p =
        LTFAT_NAME(upconv_fft_init)(L, W, a);

    LTFAT_NAME(upconv_fft_execute)(p, cin, G, F);

    LTFAT_NAME(upconv_fft_done)(p);
}

LTFAT_EXTERN LTFAT_NAME(upconv_fft_plan)
LTFAT_NAME(upconv_fft_init)(const ltfatInt L, const ltfatInt W, const ltfatInt a)
{
    ltfatInt N = L / a;
    int Nint = (int) N;

    const LTFAT_FFTW(iodim) dims = {.n = Nint, .is = 1, .os = 1};
    const LTFAT_FFTW(iodim) howmany_dims = {.n = W, .is = Nint, .os = Nint};

    LTFAT_COMPLEX* buf = ltfat_malloc(W * N * sizeof * buf);
    LTFAT_FFTW(plan) p_many =
        LTFAT_FFTW(plan_guru_dft)(1, &dims, 1, &howmany_dims,
                                  buf, buf,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

    struct LTFAT_NAME(upconv_fft_plan_struct) p_struct =
    { .L = L, .a = a, .W = W, .p_c = p_many, .buf = buf, .bufLen = W * N };

    LTFAT_NAME(upconv_fft_plan) p = ltfat_malloc(sizeof * p);
    memcpy(p, &p_struct, sizeof * p);
    return p;
}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_execute)(LTFAT_NAME(upconv_fft_plan) p,
                               const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                               LTFAT_COMPLEX *F)
{
    const ltfatInt L = p->L;
    const ltfatInt a = p->a;
    const ltfatInt W = p->W;
    LTFAT_COMPLEX* buf = p->buf;
    ltfatInt N = L / a;
    memcpy(buf, cin, W * N * sizeof * cin);


    // New array execution, inplace
    LTFAT_FFTW(execute_dft)(p->p_c, buf, buf);

    for (ltfatInt w = 0; w < W; w++)
    {
        LTFAT_COMPLEX *FPtr = F + w * L;
        LTFAT_COMPLEX *GPtr = (LTFAT_COMPLEX *) G;
        for (ltfatInt jj = 0; jj < a; jj++)
        {
            for (ltfatInt ii = 0; ii < N; ii++)
            {
                // Really readable ;)
                *FPtr++ += LTFAT_COMPLEXH(conj)(*GPtr++) * buf[ii + N * w];
            }
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_done)(LTFAT_NAME(upconv_fft_plan) p)
{
    LTFAT_FFTW(destroy_plan)(p->p_c);
    ltfat_free(p->buf);
}


LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fftbl)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                              const ltfatInt L, const ltfatInt Gl[],
                              const ltfatInt W, const double a[], const ltfatInt M,
                              const ltfatInt foff[], const int realonly[],
                              LTFAT_COMPLEX *F)
{
    // This is necessary since F us used as an accumulator
    memset(F, 0, W * L * sizeof * F);

    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_NAME(upconv_fftbl)(cin[m], G[m], L, Gl[m], W, a[m], foff[m],
                                 realonly[m], F);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fftbl_execute)(LTFAT_NAME(upconv_fftbl_plan) p[],
                                      const LTFAT_COMPLEX *cin[],
                                      const LTFAT_COMPLEX *G[],
                                      const ltfatInt M, const ltfatInt foff[],
                                      const int realonly[],
                                      LTFAT_COMPLEX *F)
{
    ltfatInt L = p[0]->L;
    ltfatInt W = p[0]->W;
    // This is necessary since F us used as an accumulator
    memset(F, 0, W * L * sizeof * F);

    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_NAME(upconv_fftbl_execute)(p[m], cin[m], G[m], foff[m], realonly[m], F);
    }

}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                         const ltfatInt L, const ltfatInt Gl, const ltfatInt W,
                         const double a,
                         const ltfatInt foff, const int realonly,
                         LTFAT_COMPLEX *F)
{
    LTFAT_NAME(upconv_fftbl_plan) p =
        LTFAT_NAME(upconv_fftbl_init)( L, Gl, W, a);

    LTFAT_NAME(upconv_fftbl_execute)(p, cin, G, foff, realonly, F);

    LTFAT_NAME(upconv_fftbl_done)( p);
}

LTFAT_EXTERN LTFAT_NAME(upconv_fftbl_plan)
LTFAT_NAME(upconv_fftbl_init)( const ltfatInt L, const ltfatInt Gl,
                               const ltfatInt W, const double a)
{
    ltfatInt N = (ltfatInt) floor(L / a + 0.5);
    int Nint = (int) N;

    const LTFAT_FFTW(iodim) dims = {.n = Nint, .is = 1, .os = 1};
    const LTFAT_FFTW(iodim) howmany_dims = {.n = W, .is = Nint, .os = Nint};

    LTFAT_COMPLEX* buf = ltfat_malloc(N * W * sizeof * buf);
    LTFAT_FFTW(plan) p_many =
        LTFAT_FFTW(plan_guru_dft)(1, &dims, 1, &howmany_dims,
                                  buf, buf,
                                  FFTW_FORWARD, FFTW_ESTIMATE);


    struct LTFAT_NAME(upconv_fftbl_plan_struct) p_struct =
    {
        .L = L, .Gl = Gl, .a = a, .W = W,
        .p_c = p_many, .buf = buf, .bufLen = N * W
    };

    LTFAT_NAME(upconv_fftbl_plan) p = ltfat_malloc(sizeof * p);
    memcpy(p, &p_struct, sizeof * p);
    return p;
}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_execute)(const LTFAT_NAME(upconv_fftbl_plan) p,
                                 const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                                 const ltfatInt foff,
                                 const int realonly, LTFAT_COMPLEX *F)
{
    const ltfatInt Gl = p->Gl;
    if(!Gl) return; // Bail out if filter has zero bandwidth
    const ltfatInt L = p->L;
    const ltfatInt W = p->W;
    const double a = p->a;
    LTFAT_COMPLEX* cbuf = p->buf;

    ltfatInt N = (ltfatInt) floor(L / a + 0.5);
    memcpy(cbuf, cin, N * W * sizeof * cin);
    LTFAT_FFTW(execute_dft)(p->p_c, cbuf, cbuf);

    for (ltfatInt w = 0; w < W; w++)
    {

        LTFAT_NAME_COMPLEX(circshift)(cbuf + w * N, cbuf + w * N, N, -foff);
        //LTFAT_NAME_COMPLEX(circshift)(cbuf+w*N,cbuf+w*N,N,Gl/2);

        const LTFAT_COMPLEX* GPtrTmp = G;
        LTFAT_COMPLEX* FPtrTmp = F + w * L;
        LTFAT_COMPLEX* CPtrTmp = cbuf + w * N;

        // Determine range of G
        ltfatInt foffTmp = foff;
        ltfatInt tmpLg = N < Gl ? N : Gl;
        ltfatInt over = 0;
        if (foffTmp + tmpLg > (ltfatInt)L)
        {
            over = foffTmp + tmpLg - (ltfatInt)L;
        }


        if (foffTmp < 0)
        {
            ltfatInt toCopy = (-foffTmp) < tmpLg ? -foffTmp : tmpLg;
            FPtrTmp = F + (w + 1) * L + foffTmp;
            for (ltfatInt ii = 0; ii < toCopy; ii++)
            {
                LTFAT_COMPLEX tmp = *CPtrTmp++ * LTFAT_COMPLEXH(conj)(*GPtrTmp++);
                FPtrTmp[ii] += tmp;
            }

            tmpLg -= toCopy;
            foffTmp = 0;
        }

        FPtrTmp = F + w * L + foffTmp;
        for (ltfatInt ii = 0; ii < tmpLg - over; ii++)
        {
            LTFAT_COMPLEX tmp = *CPtrTmp++ * LTFAT_COMPLEXH(conj)(*GPtrTmp++);
            FPtrTmp[ii] += tmp;
        }

        FPtrTmp = F + w * L;
        for (ltfatInt ii = 0; ii < over; ii++)
        {
            LTFAT_COMPLEX tmp = (*CPtrTmp++ * LTFAT_COMPLEXH(conj)(*GPtrTmp++));
            FPtrTmp[ii] += tmp;
        }
    }


    if (realonly)
    {
        const ltfatInt foffconj = -L + positiverem(L - foff - Gl, L) + 1;
        LTFAT_COMPLEX *Gconj = ltfat_malloc(Gl * sizeof * Gconj);
        LTFAT_NAME_COMPLEX(reverse_array)((LTFAT_COMPLEX *)G, Gconj, Gl);
        LTFAT_NAME_COMPLEX(conjugate_array)(Gconj, Gconj, Gl);

        LTFAT_NAME(upconv_fftbl_execute)(p, cin, Gconj, foffconj, 0, F);
        ltfat_free(Gconj);
    }

}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_done)(LTFAT_NAME(upconv_fftbl_plan) p)
{
    LTFAT_FFTW(destroy_plan)(p->p_c);
    if(p->buf) ltfat_free(p->buf);
}
