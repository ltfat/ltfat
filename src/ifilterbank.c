#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME(upconv_fft_plan_struct)
{
    ltfatInt L;
    ltfatInt W;
    ltfatInt a;
    LTFAT_FFTW(plan) p_c;
    LTFAT_COMPLEX* buf;
    ltfatInt bufLen;
};

struct LTFAT_NAME(upconv_fftbl_plan_struct)
{
    ltfatInt L;
    ltfatInt Gl;
    ltfatInt W;
    double a;
    LTFAT_FFTW(plan) p_c;
    LTFAT_COMPLEX* buf;
    ltfatInt bufLen;
};

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft)(const LTFAT_COMPLEX* cin[],
                            const LTFAT_COMPLEX* G[],
                            const ltfatInt L, const ltfatInt W, const ltfatInt a[],
                            const ltfatInt M, LTFAT_COMPLEX* F)
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
                                    const LTFAT_COMPLEX* cin[],
                                    const LTFAT_COMPLEX* G[],
                                    const ltfatInt M,
                                    LTFAT_COMPLEX* F )
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
LTFAT_NAME(upconv_fft)(const LTFAT_COMPLEX* cin, const LTFAT_COMPLEX* G,
                       const ltfatInt L, const ltfatInt W, const ltfatInt a,
                       LTFAT_COMPLEX* F)
{
    LTFAT_NAME(upconv_fft_plan) p =
        LTFAT_NAME(upconv_fft_init)(L, W, a);

    LTFAT_NAME(upconv_fft_execute)(p, cin, G, F);

    LTFAT_NAME(upconv_fft_done)(p);
}

LTFAT_EXTERN LTFAT_NAME(upconv_fft_plan)
LTFAT_NAME(upconv_fft_init)(const ltfatInt L, const ltfatInt W,
                            const ltfatInt a)
{
    ltfatInt N = L / a;

    LTFAT_FFTW(iodim64) dims;
    dims.n = N; dims.is = 1; dims.os = 1;
    LTFAT_FFTW(iodim64) howmany_dims;
    howmany_dims.n = W; howmany_dims.is = N; howmany_dims.os = N;

    LTFAT_COMPLEX* buf = LTFAT_NAME_COMPLEX(malloc)(W * N);
    LTFAT_FFTW(plan) p_many =
        LTFAT_FFTW(plan_guru64_dft)(1, &dims, 1, &howmany_dims,
                                  (LTFAT_FFTW(complex)*) buf,
                                  (LTFAT_FFTW(complex)*) buf,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

    /* struct LTFAT_NAME(upconv_fft_plan_struct) p_struct = */
    /* { .L = L, .a = a, .W = W, .p_c = p_many, .buf = buf, .bufLen = W * N }; */

    LTFAT_NAME(upconv_fft_plan) p =
        (LTFAT_NAME(upconv_fft_plan))ltfat_malloc( sizeof * p);
    p->L = L; p->a = a; p->W = W; p->p_c = p_many; p->buf = buf; p->bufLen = W * N;
    /* memcpy(p, &p_struct, sizeof * p); */
    return p;
}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_execute)(LTFAT_NAME(upconv_fft_plan) p,
                               const LTFAT_COMPLEX* cin, const LTFAT_COMPLEX* G,
                               LTFAT_COMPLEX* F)
{
    const ltfatInt L = p->L;
    const ltfatInt a = p->a;
    const ltfatInt W = p->W;
    LTFAT_COMPLEX* buf = p->buf;
    ltfatInt N = L / a;
    memcpy(buf, cin, W * N * sizeof * cin);


    // New array execution, inplace
    LTFAT_FFTW(execute_dft)(p->p_c,
                            (LTFAT_FFTW(complex)*) buf,
                            (LTFAT_FFTW(complex)*) buf);

    for (ltfatInt w = 0; w < W; w++)
    {
        LTFAT_COMPLEX* FPtr = F + w * L;
        LTFAT_COMPLEX* GPtr = (LTFAT_COMPLEX*) G;
        for (ltfatInt jj = 0; jj < a; jj++)
        {
            for (ltfatInt ii = 0; ii < N; ii++)
            {
                // Really readable ;)
                *FPtr++ += conj(*GPtr++) * buf[ii + N * w];
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
LTFAT_NAME(ifilterbank_fftbl)(const LTFAT_COMPLEX* cin[],
                              const LTFAT_COMPLEX* G[],
                              const ltfatInt L, const ltfatInt Gl[],
                              const ltfatInt W, const double a[], const ltfatInt M,
                              const ltfatInt foff[], const int realonly[],
                              LTFAT_COMPLEX* F)
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
                                      const LTFAT_COMPLEX* cin[],
                                      const LTFAT_COMPLEX* G[],
                                      const ltfatInt M, const ltfatInt foff[],
                                      const int realonly[],
                                      LTFAT_COMPLEX* F)
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
LTFAT_NAME(upconv_fftbl)(const LTFAT_COMPLEX* cin, const LTFAT_COMPLEX* G,
                         const ltfatInt L, const ltfatInt Gl, const ltfatInt W,
                         const double a,
                         const ltfatInt foff, const int realonly,
                         LTFAT_COMPLEX* F)
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
    ltfatInt bufLen =  N > Gl ? N : Gl ;

    LTFAT_FFTW(iodim64) dims;
    dims.n = N; dims.is = 1; dims.os = 1;
    LTFAT_FFTW(iodim64) howmany_dims;
    howmany_dims.n = W; howmany_dims.is = bufLen; howmany_dims.os = bufLen;

    LTFAT_COMPLEX* buf = LTFAT_NAME_COMPLEX(malloc)(bufLen * W);
    LTFAT_FFTW(plan) p_many =
        LTFAT_FFTW(plan_guru64_dft)(1, &dims, 1, &howmany_dims,
                                  (LTFAT_FFTW(complex)*)buf,
                                  (LTFAT_FFTW(complex)*)buf,
                                  FFTW_FORWARD, FFTW_ESTIMATE);


    /* struct LTFAT_NAME(upconv_fftbl_plan_struct) p_struct = */
    /* { */
    /*     .L = L, .Gl = Gl, .a = a, .W = W, */
    /*     .p_c = p_many, .buf = buf, .bufLen = bufLen */
    /* }; */

    LTFAT_NAME(upconv_fftbl_plan) p =
        (LTFAT_NAME(upconv_fftbl_plan)) ltfat_malloc(sizeof * p);
    p->L = L; p->Gl = Gl; p->a = a; p->W = W;
    p->p_c = p_many; p->buf = buf; p->bufLen = bufLen;

    /* memcpy(p, &p_struct, sizeof * p); */
    return p;
}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_execute)(const LTFAT_NAME(upconv_fftbl_plan) p,
                                 const LTFAT_COMPLEX* cin, const LTFAT_COMPLEX* G,
                                 const ltfatInt foff,
                                 const int realonly, LTFAT_COMPLEX* F)
{
    ltfatInt Gl = p->Gl;
    if (!Gl) return; // Bail out if filter has zero bandwidth
    const ltfatInt bufLen = p->bufLen;
    const ltfatInt L = p->L;
    const ltfatInt W = p->W;
    const double a = p->a;
    LTFAT_COMPLEX* cbuf = p->buf;

    ltfatInt N = (ltfatInt) floor(L / a + 0.5);

    for (ltfatInt w = 0; w < W; w++)
        memcpy(cbuf + w * bufLen, cin + w * N, N * sizeof * cin);

    LTFAT_FFTW(execute_dft)(p->p_c, (LTFAT_FFTW(complex)*)cbuf,
                            (LTFAT_FFTW(complex)*)cbuf);

    for (ltfatInt w = 0; w < W; w++)
    {

        LTFAT_NAME_COMPLEX(circshift)(cbuf + w * bufLen, N, -foff, cbuf + w * bufLen);
        // This does nothing if bufLen == N
        LTFAT_NAME_COMPLEX(periodize_array)(cbuf + w * bufLen, N, bufLen,
                                            cbuf + w * bufLen);

        const LTFAT_COMPLEX* GPtrTmp = G;
        LTFAT_COMPLEX* FPtrTmp = F + w * L;
        LTFAT_COMPLEX* CPtrTmp = cbuf + w * bufLen;
        ltfatInt Gltmp = Gl;

        // Determine range of G
        ltfatInt foffTmp = foff;

        ltfatInt over = 0;
        if (foffTmp + Gltmp > (ltfatInt)L)
        {
            over = foffTmp + Gltmp - (ltfatInt)L;
        }


        if (foffTmp < 0)
        {
            ltfatInt toCopy = (-foffTmp) < Gltmp ? -foffTmp : Gltmp;
            FPtrTmp = F + (w + 1) * L + foffTmp;
            for (ltfatInt ii = 0; ii < toCopy; ii++)
            {
                LTFAT_COMPLEX tmp = *CPtrTmp++ * conj(*GPtrTmp++);
                FPtrTmp[ii] += tmp;
            }

            Gltmp -= toCopy;
            foffTmp = 0;
        }

        FPtrTmp = F + w * L + foffTmp;
        for (ltfatInt ii = 0; ii < Gltmp - over; ii++)
        {
            LTFAT_COMPLEX tmp = *CPtrTmp++ * conj(*GPtrTmp++);
            FPtrTmp[ii] += tmp;
        }

        FPtrTmp = F + w * L;
        for (ltfatInt ii = 0; ii < over; ii++)
        {
            LTFAT_COMPLEX tmp = (*CPtrTmp++ * conj(*GPtrTmp++));
            FPtrTmp[ii] += tmp;
        }
    }


    if (realonly)
    {
        const ltfatInt foffconj = -L + ltfat_positiverem(L - foff - Gl, L) + 1;
        LTFAT_COMPLEX* Gconj = LTFAT_NAME_COMPLEX(malloc)(Gl);
        LTFAT_NAME_COMPLEX(reverse_array)(G, Gl, Gconj);
        LTFAT_NAME_COMPLEX(conjugate_array)(Gconj, Gl, Gconj);

        LTFAT_NAME(upconv_fftbl_execute)(p, cin, Gconj, foffconj, 0, F);
        ltfat_free(Gconj);
    }

}


LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_done)(LTFAT_NAME(upconv_fftbl_plan) p)
{
    LTFAT_FFTW(destroy_plan)(p->p_c);
    if (p->buf) ltfat_free(p->buf);
}
