/**/
#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "kissfft/kiss_fft.h"

/****** FFT ******/
struct LTFAT_NAME(fft_plan)
{
    ltfat_int L;
    ltfat_int W;
    LTFAT_COMPLEX* in;
    LTFAT_COMPLEX* out;
    LTFAT_COMPLEX* tmp;
    LTFAT_KISS(fft_cfg) kiss_plan;
};

LTFAT_API int
LTFAT_NAME(fft)(LTFAT_COMPLEX in[], ltfat_int L, ltfat_int W,
                LTFAT_COMPLEX out[])
{
    LTFAT_NAME(fft_plan)* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS( LTFAT_NAME(fft_init)(L, W, in, out, 0, &p),
                 "Init failed");

    LTFAT_NAME(fft_execute)(p);
    LTFAT_NAME(fft_done)(&p);
error:
    return status;
}

static int
LTFAT_NAME(fft_init_common)(ltfat_int L, ltfat_int W,
                            LTFAT_COMPLEX in[], LTFAT_COMPLEX out[],
                            unsigned inverse, LTFAT_NAME(fft_plan)** p)
{
    LTFAT_NAME(fft_plan)* fftwp = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_BADARG, ( in == NULL && out == NULL ) ||
          ( in != NULL && out != NULL ),
          "in and out must both be valid or NULL");

    CHECKMEM( fftwp = (LTFAT_NAME(fft_plan)*) ltfat_calloc(1, sizeof * fftwp) );
    fftwp->L = L; fftwp->W = W; fftwp->in = in; fftwp->out = out;

    fftwp->kiss_plan = LTFAT_KISS(fft_alloc)(L, inverse, NULL, NULL);
    CHECKINIT(fftwp->kiss_plan, "FFTW plan creation failed.");

    if (in == out)
        CHECKMEM( fftwp->tmp = LTFAT_NAME_COMPLEX(malloc)(L) );

    *p = fftwp;
    return status;
error:
    LTFAT_NAME(fft_done)(&fftwp);
    *p = NULL;
    return status;
}


LTFAT_API int
LTFAT_NAME(fft_init)(ltfat_int L, ltfat_int W,
                     LTFAT_COMPLEX in[], LTFAT_COMPLEX out[],
                     unsigned UNUSED(flags), LTFAT_NAME(fft_plan)** p)
{
    return LTFAT_NAME(fft_init_common)(L, W, in, out, 0, p);
}

LTFAT_API int
LTFAT_NAME(fft_execute)(LTFAT_NAME(fft_plan)* p)
{
    return LTFAT_NAME(fft_execute_newarray)( p, p->in, p->out);
}

LTFAT_API int
LTFAT_NAME(fft_execute_newarray)(LTFAT_NAME(fft_plan)* p,
                                 const LTFAT_COMPLEX in[], LTFAT_COMPLEX out[])
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(in); CHECKNULL(out);

    if (in == out)
    {
        CHECKNULL(p->tmp);

        for (ltfat_int w = 0; w < p->W; w++)
        {
            memcpy(p->tmp, in + w * p->L, p->L * sizeof * p->tmp);
            LTFAT_KISS(fft)(p->kiss_plan,
                            (const kiss_fft_cpx*) p->tmp,
                            (kiss_fft_cpx*) out + w * p->L);
        }
    }
    else
    {
        for (ltfat_int w = 0; w < p->W; w++)
            LTFAT_KISS(fft)(p->kiss_plan,
                            (const kiss_fft_cpx*) in + w * p->L,
                            (kiss_fft_cpx*) out + w * p->L);
    }

error:
    return status;
}

LTFAT_API int
LTFAT_NAME(fft_done)(LTFAT_NAME(fft_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    LTFAT_NAME(fft_plan)* pp = NULL;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    if (pp->tmp) ltfat_free(pp->tmp);
    if (pp->kiss_plan) ltfat_free(pp->kiss_plan);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

/******* IFFT ******/
struct LTFAT_NAME(ifft_plan)
{
    struct LTFAT_NAME(fft_plan) inplan;
};

LTFAT_API int
LTFAT_NAME(ifft)(LTFAT_COMPLEX in[], ltfat_int L, ltfat_int W,
                 LTFAT_COMPLEX out[])
{
    LTFAT_NAME(ifft_plan)* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS( LTFAT_NAME(ifft_init)(L, W, in, out, 0, &p),
                 "Init failed");

    LTFAT_NAME(ifft_execute)(p);
    LTFAT_NAME(ifft_done)(&p);
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(ifft_init)(ltfat_int L, ltfat_int W,
                      LTFAT_COMPLEX in[], LTFAT_COMPLEX out[],
                      unsigned UNUSED(flags), LTFAT_NAME(ifft_plan)** p)
{
    return LTFAT_NAME(fft_init_common)(L, W, in, out, 1,
                                       (LTFAT_NAME(fft_plan)**) p);
}

LTFAT_API int
LTFAT_NAME(ifft_execute)(LTFAT_NAME(ifft_plan)* p)
{
    return LTFAT_NAME(fft_execute)((LTFAT_NAME(fft_plan)*) p);
}

LTFAT_API int
LTFAT_NAME(ifft_execute_newarray)(LTFAT_NAME(ifft_plan)* p,
                                  const LTFAT_COMPLEX in[], LTFAT_COMPLEX out[])
{
    return LTFAT_NAME(fft_execute_newarray)((LTFAT_NAME(fft_plan)*) p, in, out);

}

LTFAT_API int
LTFAT_NAME(ifft_done)(LTFAT_NAME(ifft_plan)** p)
{
    return LTFAT_NAME(fft_done)((LTFAT_NAME(fft_plan)**) p);
}

/****** FFTREAL ******/
struct LTFAT_NAME(fftreal_plan)
{
    ltfat_int L;
    ltfat_int W;
    LTFAT_REAL* in;
    LTFAT_REAL* out;
    LTFAT_COMPLEX* tmp;
    LTFAT_KISS(fft_cfg) kiss_plan;
    LTFAT_KISS(fftr_cfg) kiss_plan;
};

static int
LTFAT_NAME(fftreal_init_common)(ltfat_int L, ltfat_int W,
                                LTFAT_REAL in[], LTFAT_REAL out[],
                                unsigned inverse, LTFAT_NAME(fftreal_plan)** p)
{
    LTFAT_NAME(fftreal_plan)* fftwp = NULL;

    int status = LTFATERR_SUCCESS;

    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, L > 0 && ~(L & 1), "L must be even positive");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");
    CHECK(LTFATERR_BADARG, ( in == NULL && out == NULL ) ||
          ( in != NULL && out != NULL ),
          "in and out must both be valid or NULL");

    CHECKMEM( fftwp = (LTFAT_NAME(fftreal_plan)*) ltfat_calloc(1, sizeof * fftwp) );
    fftwp->L = L; fftwp->W = W; fftwp->in = in; fftwp->out = out;

    fftwp->kiss_plan = LTFAT_KISS(fftr_alloc)(L, inverse, NULL, NULL);
    CHECKINIT(fftwp->kiss_plan, "FFTW plan creation failed.");

    if (in == out)
        CHECKMEM( fftwp->tmp = LTFAT_NAME_COMPLEX(malloc)(L / 2 + 1) );

    *p = fftwp;
    return status;
error:
    LTFAT_NAME(fftreal_done)(&fftwp);
    *p = NULL;
    return status;
}


LTFAT_API int
LTFAT_NAME(fftreal_init)(ltfat_int L, ltfat_int W,
                         LTFAT_REAL in[], LTFAT_COMPLEX out[],
                         unsigned UNUSED(flags), LTFAT_NAME(fftreal_plan)** p)
{
    if (L & 1)
        return  LTFAT_NAME(fft_init_common)(L, W, in, (LTFAT_REAL*) out, 0,
                                            (LTFAT_NAME(fft_plan)**) p);
    else
        return  LTFAT_NAME(fftreal_init_common)(L, W, in, (LTFAT_REAL*) out, 0, p);
}

LTFAT_API int
LTFAT_NAME(fftreal_execute)(LTFAT_NAME(fftreal_plan)* p)
{
    return LTFAT_NAME(fftreal_execute_newarray)( p, p->in,
            (LTFAT_COMPLEX*) p->out);
}

LTFAT_API int
LTFAT_NAME(fftreal_execute_newarray)(LTFAT_NAME(fftreal_plan)* p,
                                     const LTFAT_REAL in[], LTFAT_COMPLEX out[])
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(in); CHECKNULL(out);
    ltfat_int M2 = p->L / 2 + 1;

    if ((in == (const LTFAT_REAL*) out) || (p->L & 1) )
    {
        CHECKNULL(p->tmp);

        if (p->L & 1)
        {
            for (ltfat_int w = 0; w < p->W; w++)
            {
                LTFAT_NAME(real2complex_array)(in + w*p->L, L, p->tmp);
                LTFAT_KISS(fft)(p->kiss_plan,
                                (const kiss_fft_scalar*) p->tmp,
                                (kiss_fft_cpx*) out + w * M2);
            }
        }
        else
        {
            for (ltfat_int w = 0; w < p->W; w++)
            {
                memcpy(p->tmp, in + w * p->L, p->L * sizeof * p->in);
                LTFAT_KISS(fftr)(p->kiss_plan,
                                 (const kiss_fft_scalar*) p->tmp,
                                 (kiss_fft_cpx*) out + w * M2);
            }
        }
    }
    else
    {
        for (ltfat_int w = 0; w < p->W; w++)
            LTFAT_KISS(fftr)(p->kiss_plan,
                             (const kiss_fft_scalar*) in + w * p->L,
                             (kiss_fft_cpx*) out + w * M2);
    }
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(fftreal_done)(LTFAT_NAME(fftreal_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    LTFAT_NAME(fftreal_plan)* pp = NULL;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    if (pp->tmp) ltfat_free(pp->tmp);
    if (pp->kiss_plan) ltfat_free(pp->kiss_plan);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

/******* IFFTREAL ******/
struct LTFAT_NAME(ifftreal_plan)
{
    LTFAT_NAME(fftreal_plan) inplan;
};

LTFAT_API int
LTFAT_NAME(ifftreal)(LTFAT_COMPLEX in[], ltfat_int L, ltfat_int W,
                     LTFAT_REAL out[])
{
    LTFAT_NAME(ifftreal_plan)* p = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKSTATUS( LTFAT_NAME(ifftreal_init)(L, W, in, out, 0, &p),
                 "Init failed");

    LTFAT_NAME(ifftreal_execute)(p);
LTFAT_NAME(ifftreal_done)(&p); error:
    return status;
}

LTFAT_API int
LTFAT_NAME(ifftreal_init)(ltfat_int L, ltfat_int W,
                          LTFAT_COMPLEX in[], LTFAT_REAL out[],
                          unsigned UNUSED(flags), LTFAT_NAME(ifftreal_plan)** p)
{
    return LTFAT_NAME(fftreal_init_common)(L, W, (LTFAT_REAL*)in, out, 1,
                                           (LTFAT_NAME(fftreal_plan)**) p);
}

LTFAT_API int
LTFAT_NAME(ifftreal_execute)(LTFAT_NAME(ifftreal_plan)* pin)
{
    LTFAT_NAME(fftreal_plan)* p = (LTFAT_NAME(fftreal_plan)*) pin;
    return LTFAT_NAME(ifftreal_execute_newarray)( pin, (const LTFAT_COMPLEX*) p->in,
            p->out);
}

LTFAT_API int
LTFAT_NAME(ifftreal_execute_newarray)(LTFAT_NAME(ifftreal_plan)* pin,
                                      const LTFAT_COMPLEX in[], LTFAT_REAL out[])
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(pin); CHECKNULL(in); CHECKNULL(out);
    LTFAT_NAME(fftreal_plan)* p = (LTFAT_NAME(fftreal_plan)*) pin;

    ltfat_int M2 = p->L / 2 + 1;

    if (in == (const LTFAT_COMPLEX*) out)
    {
        CHECKNULL(p->tmp);

        for (ltfat_int w = 0; w < p->W; w++)
        {
            memcpy(p->tmp, in + w * M2, M2 * sizeof * p->in);
            LTFAT_KISS(fftri)(p->kiss_plan,
                              (const kiss_fft_cpx*) p->tmp,
                              (kiss_fft_scalar*) out + w * p->L);
        }
    }
    else
    {
        for (ltfat_int w = 0; w < p->W; w++)
            LTFAT_KISS(fftri)(p->kiss_plan,
                              (const kiss_fft_cpx*) in + w * M2,
                              (kiss_fft_scalar*) out + w * p->L);
    }
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(ifftreal_done)(LTFAT_NAME(ifftreal_plan)** p)
{
    return LTFAT_NAME(fftreal_done)((LTFAT_NAME(fftreal_plan)**) p);
}
