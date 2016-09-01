#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

// These are non-public function header templates
typedef int LTFAT_NAME(realtocomplextransform)(void* userdata,
        const LTFAT_REAL* in, ltfat_int, LTFAT_COMPLEX* out);

typedef int LTFAT_NAME(complextorealtransform)(void* userdata,
        const LTFAT_COMPLEX* in, ltfat_int W, LTFAT_REAL* out);

struct LTFAT_NAME(rtdgtreal_plan)
{
    LTFAT_REAL* g; //!< Window
    ltfat_int gl; //!< Window length
    ltfat_int M; //!< Number of FFT channels
    rtdgt_phasetype ptype; //!< Phase convention
    LTFAT_REAL* fftBuf; //!< Internal buffer
    ltfat_int fftBufLen; //!< Internal buffer length
    LTFAT_FFTW(plan) pfft; //!< FFTW plan
};


int
LTFAT_NAME(rtdgtreal_commoninit)(const LTFAT_REAL* g, ltfat_int gl,
                                 ltfat_int M, const rtdgt_phasetype ptype,
                                 const ltfat_transformdirection tradir,
                                 LTFAT_NAME(rtdgtreal_plan)** pout)
{
    ltfat_int M2;
    LTFAT_NAME(rtdgtreal_plan)* p = NULL;
    LTFAT_FFTW(iodim64) dims;
    LTFAT_FFTW(iodim64) howmany_dims;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    CHECKMEM( p = (LTFAT_NAME(rtdgtreal_plan)*) ltfat_calloc(1, sizeof * p) );

    M2 = M / 2 + 1;
    p->fftBufLen = gl > 2 * M2 ? gl : 2 * M2;

    CHECKMEM( p->g = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM( p->fftBuf = LTFAT_NAME_REAL(malloc)(p->fftBufLen));
    p->gl = gl;
    p->M = M;
    p->ptype = ptype;

    LTFAT_NAME_REAL(fftshift)(g, gl, p->g);

    if (LTFAT_FORWARD == tradir)
    {
        dims.n = M; dims.is = 1; dims.os = 1;
        howmany_dims.n = 1; howmany_dims.is = M; howmany_dims.os = M / 2 + 1;

        p->pfft =
            LTFAT_FFTW(plan_guru64_dft_r2c)(1, &dims, 1, &howmany_dims,
                                            p->fftBuf, (LTFAT_FFTW(complex)*) p->fftBuf,
                                            FFTW_MEASURE);
    }
    else if (LTFAT_INVERSE == tradir)
    {
        dims.n = M; dims.is = 1; dims.os = 1;
        howmany_dims.n = 1; howmany_dims.is = M / 2 + 1; howmany_dims.os = M;

        p->pfft =
            LTFAT_FFTW(plan_guru64_dft_c2r)(1, &dims, 1, &howmany_dims,
                                            (LTFAT_FFTW(complex)*)
                                            p->fftBuf, p->fftBuf, FFTW_MEASURE);
    }
    else
        CHECKCANTHAPPEN("Unknown transform direction.");

    CHECKINIT(p->pfft, "FFTW plan creation failed.");

    *pout = p;
    return status;
error:
    if (p) LTFAT_NAME(rtdgtreal_done)(&p);
    return status;
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_init)(const LTFAT_REAL* g, ltfat_int gl,
                           ltfat_int M, const rtdgt_phasetype ptype,
                           LTFAT_NAME(rtdgtreal_plan)** p)
{
    return LTFAT_NAME(rtdgtreal_commoninit)(g, gl, M, ptype, LTFAT_FORWARD, p);
}

LTFAT_API int
LTFAT_NAME(rtidgtreal_init)(const LTFAT_REAL* g, ltfat_int gl,
                            ltfat_int M, const rtdgt_phasetype ptype,
                            LTFAT_NAME(rtdgtreal_plan)** p)
{
    return LTFAT_NAME(rtdgtreal_commoninit)(g, gl, M, ptype, LTFAT_INVERSE, p);
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_execute)(const LTFAT_NAME(rtdgtreal_plan)* p,
                              const LTFAT_REAL* f, ltfat_int W,
                              LTFAT_COMPLEX* c)
{
    ltfat_int M, M2, gl;
    LTFAT_REAL* fftBuf;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(f); CHECKNULL(c);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");

    M = p->M;
    M2 = M / 2 + 1;
    gl = p->gl;
    fftBuf = p->fftBuf;

    for (ltfat_int w = 0; w < W; w++)
    {
        const LTFAT_REAL* fchan = f + w * gl;
        LTFAT_COMPLEX* cchan = c + w * M2;

        if (p->g)
            for (ltfat_int ii = 0; ii < gl; ii++)
                fftBuf[ii] = fchan[ii] * p->g[ii];

        if (M > gl)
            memset(fftBuf + gl, 0, (M - gl) * sizeof * fftBuf);

        if (gl > M)
            LTFAT_NAME_REAL(fold_array)(fftBuf, gl, M, 0, fftBuf);

        if (p->ptype == LTFAT_RTDGTPHASE_ZERO)
            LTFAT_NAME_REAL(circshift)(fftBuf, M, -(gl / 2), fftBuf );

        LTFAT_FFTW(execute)(p->pfft);

        memcpy(cchan, fftBuf, M2 * sizeof * c);
    }

error:
    return status;
}

int
LTFAT_NAME(rtdgtreal_execute_wrapper)(void* p,
                                      const LTFAT_REAL* f, ltfat_int W,
                                      LTFAT_COMPLEX* c)
{
    return LTFAT_NAME(rtdgtreal_execute)((LTFAT_NAME(rtdgtreal_plan)*) p, f, W, c);
}

LTFAT_API int
LTFAT_NAME(rtidgtreal_execute)(const LTFAT_NAME(rtidgtreal_plan)* p,
                               const LTFAT_COMPLEX* c, ltfat_int W,
                               LTFAT_REAL* f)
{
    ltfat_int M, M2, gl;
    LTFAT_REAL* fftBuf;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(c); CHECKNULL(f);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive");

    M = p->M;
    M2 = M / 2 + 1;
    gl = p->gl;
    fftBuf = p->fftBuf;

    for (ltfat_int w = 0; w < W; w++)
    {
        const LTFAT_COMPLEX* cchan = c + w * M2;
        LTFAT_REAL* fchan = f + w * gl;

        memcpy(fftBuf, cchan, M2 * sizeof * cchan);

        LTFAT_FFTW(execute)(p->pfft);

        if (p->ptype == LTFAT_RTDGTPHASE_ZERO)
            LTFAT_NAME_REAL(circshift)(fftBuf, M, gl / 2, fftBuf );

        if (gl > M)
            LTFAT_NAME_REAL(periodize_array)(fftBuf, M , gl, fftBuf);

        if (p->g)
            for (ltfat_int ii = 0; ii < gl; ii++)
                fftBuf[ii] *= p->g[ii];

        memcpy(fchan, fftBuf, gl * sizeof * fchan);
    }
error:
    return status;
}

int
LTFAT_NAME(rtidgtreal_execute_wrapper)(void* p,
                                       const LTFAT_COMPLEX* c, ltfat_int W,
                                       LTFAT_REAL* f)
{
    return LTFAT_NAME(rtidgtreal_execute)((LTFAT_NAME(rtidgtreal_plan)*)p, c, W, f);
}


LTFAT_API int
LTFAT_NAME(rtdgtreal_done)(LTFAT_NAME(rtdgtreal_plan)** p)
{
    LTFAT_NAME(rtdgtreal_plan)* pp;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);

    pp = *p;
    ltfat_safefree(pp->g);
    ltfat_safefree(pp->fftBuf);
    if (pp->pfft) LTFAT_FFTW(destroy_plan)(pp->pfft);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(rtidgtreal_done)(LTFAT_NAME(rtidgtreal_plan)** p)
{
    return LTFAT_NAME(rtdgtreal_done)(p);
}

/* FWD FIFO */

struct LTFAT_NAME(rtdgtreal_fifo_state)
{
    ltfat_int gl; //!< Window length
    ltfat_int a; //!< Hop size
    LTFAT_REAL* buf; //!< Ring buffer array
    ltfat_int bufLen; //!< Length of the previous
    ltfat_int readIdx; //!< Read pos.
    ltfat_int writeIdx; //!< Write pos.
    ltfat_int Wmax;
};


LTFAT_API int
LTFAT_NAME(rtdgtreal_fifo_init)(ltfat_int fifoLen, ltfat_int procDelay,
                                ltfat_int gl, ltfat_int a, ltfat_int Wmax,
                                LTFAT_NAME(rtdgtreal_fifo_state)** pout)
{
    LTFAT_NAME(rtdgtreal_fifo_state)* p = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, fifoLen > 0, "fifoLen must be positive");
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");
    CHECK(LTFATERR_BADARG, procDelay >= gl - 1 , "procDelay must be positive");
    CHECK(LTFATERR_BADARG , fifoLen > gl + 1, "fifoLen must be bugger than gl+1");

    CHECKMEM(p = (LTFAT_NAME(rtdgtreal_fifo_state)*) ltfat_calloc(1, sizeof * p));
    CHECKMEM(p->buf = LTFAT_NAME_REAL(calloc)( Wmax * (fifoLen + 1)));

    p->bufLen = fifoLen + 1;
    p->a = a; p->gl = gl; p->readIdx = fifoLen + 1 - (procDelay); p->Wmax = Wmax;

    *pout = p;
    return status;
error:
    if (p) LTFAT_NAME(rtdgtreal_fifo_done)(&p);
    return status;
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_fifo_done)(LTFAT_NAME(rtdgtreal_fifo_state)** p)
{
    LTFAT_NAME(rtdgtreal_fifo_state)* pp;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    ltfat_safefree(pp->buf);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_fifo_reset)(LTFAT_NAME(rtdgtreal_fifo_state)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    memset(p->buf, 0, p->Wmax * p->bufLen * sizeof * p->buf);

error:
    return status;
}

LTFAT_API ltfat_int
LTFAT_NAME(rtdgtreal_fifo_write)(LTFAT_NAME(rtdgtreal_fifo_state)* p,
                                 const LTFAT_REAL** buf,
                                 ltfat_int bufLen, ltfat_int W)
{
    ltfat_int Wact, freeSpace, toWrite, valid, over, endWriteIdx;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(buf);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive.");
    for (ltfat_int w = 0; w < W; w++)
        CHECKNULL(buf[w]);

    CHECK(LTFATERR_NOTPOSARG, bufLen > 0, "bufLen must be positive.");

    freeSpace = p->readIdx - p->writeIdx - 1;
    if (freeSpace < 0) freeSpace += p->bufLen;

    // CHECK(LTFATERR_OVERFLOW, freeSpace, "FIFO owerflow");

    Wact = p->Wmax < W ? p->Wmax : W;

    toWrite = bufLen > freeSpace ? freeSpace : bufLen;
    valid = toWrite;
    over = 0;

    endWriteIdx = p->writeIdx + toWrite;

    if (endWriteIdx > p->bufLen)
    {
        valid = p->bufLen - p->writeIdx;
        over = endWriteIdx - p->bufLen;
    }

    if (valid > 0)
    {
        for (ltfat_int w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen + p->writeIdx;
            if (w < Wact)
                memcpy(pbufchan, buf[w], valid * sizeof * p->buf );
            else
                memset(pbufchan, 0, valid * sizeof * p->buf );
        }
    }
    if (over > 0)
    {
        for (ltfat_int w = 0; w < Wact; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen;
            if (w < Wact)
                memcpy(pbufchan, buf[w] + valid, over * sizeof * p->buf);
            else
                memset(pbufchan, 0,  over * sizeof * p->buf);
        }
    }
    p->writeIdx = ( p->writeIdx + toWrite ) % p->bufLen;

    return toWrite;
error:
    return status;
}

LTFAT_API ltfat_int
LTFAT_NAME(rtdgtreal_fifo_read)(LTFAT_NAME(rtdgtreal_fifo_state)* p,
                                LTFAT_REAL* buf)
{
    ltfat_int available, toRead, valid, over, endReadIdx;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(buf);

    available = p->writeIdx - p->readIdx;
    if (available < 0) available += p->bufLen;

    // CHECK(LTFATERR_UNDERFLOW, available >= p->gl, "FIFO underflow");
    if (available < p->gl) return 0;

    toRead = p->gl;

    valid = toRead;
    over = 0;

    endReadIdx = p->readIdx + valid;

    if (endReadIdx > p->bufLen)
    {
        valid = p->bufLen - p->readIdx;
        over = endReadIdx - p->bufLen;
    }

    if (valid > 0)
    {
        for (ltfat_int w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen + p->readIdx;
            memcpy(buf + w * p->gl, pbufchan, valid * sizeof * p->buf );
        }
    }
    if (over > 0)
    {
        for (ltfat_int w = 0; w < p->Wmax; w++)
        {
            memcpy(buf + valid + w * p->gl, p->buf + w * p->bufLen, over * sizeof * p->buf);
        }
    }

    // Only advance by a
    p->readIdx = ( p->readIdx + p->a ) % p->bufLen;

    return toRead;
error:
    return status;
}

/* BACK FIFO */

struct LTFAT_NAME(rtidgtreal_fifo_state)
{
    ltfat_int gl; //!< Window length
    ltfat_int a; //!< Hop size
    LTFAT_REAL* buf; //!< Ring buffer array
    ltfat_int bufLen; //!< Length of the previous
    ltfat_int readIdx; //!< Read pos.
    ltfat_int writeIdx; //!< Write pos.
    ltfat_int Wmax;
};


LTFAT_API int
LTFAT_NAME(rtidgtreal_fifo_init)(ltfat_int fifoLen, ltfat_int gl,
                                 ltfat_int a, ltfat_int Wmax,
                                 LTFAT_NAME(rtidgtreal_fifo_state)** pout)
{
    LTFAT_NAME(rtidgtreal_fifo_state)* p = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, fifoLen > 0, "fifoLen must be positive");
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");
    CHECK(LTFATERR_BADARG , fifoLen > gl + 1, "fifoLen must be bugger than gl+1");

    CHECKMEM( p = (LTFAT_NAME(rtidgtreal_fifo_state)*) ltfat_calloc(1, sizeof * p));
    CHECKMEM( p->buf = LTFAT_NAME_REAL(calloc)( Wmax * (fifoLen + gl + 1)));
    p->a = a; p->gl = gl; p->Wmax = Wmax; p->bufLen = fifoLen + gl + 1;

    *pout = p;
    return status;
error:
    if (p) LTFAT_NAME(rtidgtreal_fifo_done)(&p);
    return status;
}

LTFAT_API int
LTFAT_NAME(rtidgtreal_fifo_reset)(LTFAT_NAME(rtidgtreal_fifo_state)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    memset(p->buf, 0, p->Wmax * p->bufLen * sizeof * p->buf);

error:
    return status;
}

LTFAT_API int
LTFAT_NAME(rtidgtreal_fifo_done)(LTFAT_NAME(rtidgtreal_fifo_state)** p)
{
    LTFAT_NAME(rtidgtreal_fifo_state)* pp;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    ltfat_safefree(pp->buf);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

LTFAT_API ltfat_int
LTFAT_NAME(rtidgtreal_fifo_write)(LTFAT_NAME(rtidgtreal_fifo_state)* p,
                                  const LTFAT_REAL* buf)
{
    ltfat_int freeSpace, toWrite, valid, over, endWriteIdx;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(buf);

    freeSpace = p->readIdx - p->writeIdx - 1;
    if (freeSpace < 0) freeSpace += p->bufLen;

    // CHECK(LTFATERR_OVERFLOW, freeSpace >= p->gl, "FIFO overflow");
    if (freeSpace < p->gl) return 0;

    toWrite = p->gl;
    valid = toWrite;
    over = 0;

    endWriteIdx = p->writeIdx + toWrite;

    if (endWriteIdx > p->bufLen)
    {
        valid = p->bufLen - p->writeIdx;
        over = endWriteIdx - p->bufLen;
    }

    if (valid > 0)
    {
        for (ltfat_int w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + p->writeIdx + w * p->bufLen;
            const LTFAT_REAL* bufchan = buf + w * p->gl;
            for (ltfat_int ii = 0; ii < valid; ii++)
                pbufchan[ii] += bufchan[ii];
        }
    }
    if (over > 0)
    {
        for (ltfat_int w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen;
            const LTFAT_REAL* bufchan = buf + valid + w * p->gl;
            for (ltfat_int ii = 0; ii < over; ii++)
                pbufchan[ii] += bufchan[ii];
        }
    }

    p->writeIdx = ( p->writeIdx + p->a ) % p->bufLen;

    return toWrite;
error:
    return status;
}

LTFAT_API ltfat_int
LTFAT_NAME(rtidgtreal_fifo_read)(LTFAT_NAME(rtidgtreal_fifo_state)* p,
                                 ltfat_int bufLen, ltfat_int W,
                                 LTFAT_REAL** buf)
{
    ltfat_int available, toRead, valid, over, endReadIdx;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(buf);
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive.");

    for (ltfat_int w = 0; w < W; w++)
        CHECKNULL(buf[w]);

    CHECK(LTFATERR_NOTPOSARG, bufLen > 0, "bufLen must be positive.");

    available = p->writeIdx - p->readIdx;
    if (available < 0) available += p->bufLen;

    // CHECK(LTFATERR_UNDERFLOW, available, "FIFO underflow");

    toRead = available < bufLen ? available : bufLen;

    valid = toRead;
    over = 0;

    endReadIdx = p->readIdx + valid;

    if (endReadIdx > p->bufLen)
    {
        valid = p->bufLen - p->readIdx;
        over = endReadIdx - p->bufLen;
    }

    // Set the just read samples to zero so that the values are not used in
    // write again
    if (valid > 0)
    {
        for (ltfat_int w = 0; w < W; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + p->readIdx + w * p->bufLen;
            memcpy(buf[w], pbufchan, valid * sizeof * p->buf);
            memset(pbufchan, 0, valid * sizeof * p->buf);
        }
    }
    if (over > 0)
    {
        for (ltfat_int w = 0; w < W; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen;
            memcpy(buf[w] + valid, pbufchan, over * sizeof * p->buf);
            memset(pbufchan, 0, over * sizeof * p->buf);
        }
    }

    p->readIdx = ( p->readIdx + toRead ) % p->bufLen;

    return toRead;
error:
    return status;
}

/* DGTREAL processor */
struct LTFAT_NAME(rtdgtreal_processor_state)
{
    LTFAT_NAME(rtdgtreal_processor_callback)*
    processorCallback; //!< Custom processor callback
    void* userdata; //!< Callback data
    LTFAT_NAME(rtdgtreal_fifo_state)* fwdfifo;
    LTFAT_NAME(rtidgtreal_fifo_state)* backfifo;
    LTFAT_NAME(rtdgtreal_plan)* fwdplan;
    LTFAT_NAME(rtidgtreal_plan)* backplan;
    LTFAT_NAME(realtocomplextransform)* fwdtra;
    LTFAT_NAME(complextorealtransform)* backtra;
    LTFAT_REAL* buf;
    LTFAT_COMPLEX* fftbufIn;
    LTFAT_COMPLEX* fftbufOut;
    ltfat_int bufLenMax;
    void** garbageBin;
    int garbageBinSize;
    const LTFAT_REAL** inTmp;
    LTFAT_REAL** outTmp;
};


LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_init)(const LTFAT_REAL* ga, ltfat_int gal,
                                     const LTFAT_REAL* gs, ltfat_int gsl,
                                     ltfat_int a, ltfat_int M, ltfat_int Wmax,
                                     ltfat_int bufLenMax,
                                     LTFAT_NAME(rtdgtreal_processor_state)** pout)
{
    LTFAT_NAME(rtdgtreal_processor_state)* p = NULL;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(pout);
    CHECK(LTFATERR_BADSIZE, gal > 0, "gla must be positive");
    CHECK(LTFATERR_BADSIZE, gsl > 0, "gls must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");
    CHECK(LTFATERR_NOTPOSARG, bufLenMax > 0, "bufLenMax must be positive");
    CHECKMEM( p =
                  (LTFAT_NAME(rtdgtreal_processor_state)*) ltfat_calloc(1, sizeof * p));

    CHECKMEM(
        p->fftbufIn = LTFAT_NAME_COMPLEX(malloc)( Wmax * (M / 2 + 1)));
    CHECKMEM(
        p->fftbufOut = LTFAT_NAME_COMPLEX(malloc)( Wmax * (M / 2 + 1)));
    CHECKMEM( p->buf = LTFAT_NAME_REAL(malloc)( Wmax * gal));

    CHECKMEM( p->inTmp = (const LTFAT_REAL**) ltfat_malloc(Wmax * sizeof *
                         p->inTmp));
    CHECKMEM( p->outTmp = (LTFAT_REAL**) ltfat_malloc(Wmax * sizeof * p->outTmp));

    CHECKSTATUS(
        LTFAT_NAME(rtdgtreal_fifo_init)(bufLenMax + gal, gal > gsl ? gal - 1 : gsl - 1 ,
                                        gal,
                                        a, Wmax, &p->fwdfifo), "fwd fifo init failed");

    CHECKSTATUS(
        LTFAT_NAME(rtidgtreal_fifo_init)(bufLenMax + gsl, gsl, a, Wmax, &p->backfifo),
        "back fifo init failed");

    CHECKSTATUS( LTFAT_NAME(rtdgtreal_init)(ga, gal, M, LTFAT_RTDGTPHASE_ZERO,
                                            &p->fwdplan), "fwd plan failed");

    CHECKSTATUS( LTFAT_NAME(rtidgtreal_init)(gs, gsl, M,
                 LTFAT_RTDGTPHASE_ZERO, &p->backplan), "back plan failed");

    p->fwdtra = &LTFAT_NAME(rtdgtreal_execute_wrapper);
    p->backtra = &LTFAT_NAME(rtidgtreal_execute_wrapper);
    p->bufLenMax = bufLenMax;

    *pout = p;
    return status;
error:
    if (p) LTFAT_NAME(rtdgtreal_processor_done)(&p);
    return status;
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_reset)(LTFAT_NAME(rtdgtreal_processor_state)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    LTFAT_NAME(rtdgtreal_fifo_reset)(p->fwdfifo);
    LTFAT_NAME(rtidgtreal_fifo_reset)(p->backfifo);

error:
    return status;
}


LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_init_win)(LTFAT_FIRWIN win,
        ltfat_int gl, ltfat_int a, ltfat_int M,
        ltfat_int Wmax, ltfat_int bufLenMax,
        LTFAT_NAME(rtdgtreal_processor_state)** pout)
{
    LTFAT_NAME(rtdgtreal_processor_state)* p;
    LTFAT_REAL* g = NULL;
    LTFAT_REAL* gd = NULL;
    void** garbageBin = NULL;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(pout);
    CHECK(LTFATERR_BADSIZE, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0,  "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0,  "M must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");

    CHECKMEM(g = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM(gd = LTFAT_NAME_REAL(malloc)(gl));
    CHECKMEM(garbageBin = (void**) ltfat_malloc(2 * sizeof(void*)));

    CHECKSTATUS(LTFAT_NAME_REAL(firwin)(win, gl, g), "Call to firwin failed");
    CHECKSTATUS(LTFAT_NAME_REAL(gabdual_painless)(g, gl, a, M, gd),
                "Call to gabdual_painless failed");

    CHECKSTATUS(LTFAT_NAME(rtdgtreal_processor_init)(g, gl, gd, gl, a, M, Wmax,
                bufLenMax, pout), "processor_init failed");

    p = *pout;
    p->garbageBinSize = 2;
    p->garbageBin = garbageBin;
    p->garbageBin[0] = g;
    p->garbageBin[1] = gd;

    return status;
error:
    LTFAT_SAFEFREEALL(g, gd, garbageBin);
    // Also status is now set to the proper value
    return status;
}


LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_setcallback)(
    LTFAT_NAME( rtdgtreal_processor_state)* p,
    LTFAT_NAME(rtdgtreal_processor_callback)* callback,
    void* userdata)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    p->processorCallback = callback;
    p->userdata = userdata;
error:
    return status;
}


LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_execute_compact)(
    LTFAT_NAME(rtdgtreal_processor_state)* p, const LTFAT_REAL* in,
    ltfat_int len, ltfat_int chanNo, LTFAT_REAL* out)
{
    ltfat_int chanLoc;
    int status2 = LTFATERR_SUCCESS;
    int status = LTFATERR_SUCCESS;

    CHECKNULL(p);
    chanLoc = chanNo > p->fwdfifo->Wmax ? p->fwdfifo->Wmax : chanNo;

    for (ltfat_int w = 0; w < chanLoc; w++)
    {
        p->inTmp[w] = &in[w * len];
        p->outTmp[w] = &out[w * len];
    }

    // Clear superfluous channels
    if (chanNo > chanLoc)
    {
        DEBUG("Channel overflow (passed %td, max %td)", chanNo, chanLoc);
        status = LTFATERR_OVERFLOW;

        memset(out + chanLoc * len, 0, (chanNo - chanLoc)*len * sizeof * out);
    }

    status2 = LTFAT_NAME(rtdgtreal_processor_execute)( p, p->inTmp, len, chanLoc,
              p->outTmp);

    if (status2 != LTFATERR_SUCCESS) return status2;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_execute)(
    LTFAT_NAME(rtdgtreal_processor_state)* p,
    const LTFAT_REAL** in, ltfat_int len, ltfat_int chanNo,
    LTFAT_REAL** out)
{
    int status = LTFATERR_SUCCESS;
    ltfat_int samplesWritten = 0, samplesRead = 0;
    // Get default processor if none was set
    LTFAT_NAME(rtdgtreal_processor_callback)* processorCallback =
        p->processorCallback;

    // Failing these checks prohibits execution altogether
    CHECKNULL(p); CHECKNULL(in); CHECKNULL(out);
    CHECK(LTFATERR_BADSIZE,
          len >= 0, "len must be positive or zero (passed %td)", len);
    CHECK(LTFATERR_BADSIZE,
          chanNo >= 0, "chanNo must be positive or zero (passed %td)", chanNo);

    // Just dont do anything
    if (len == 0 || chanNo == 0) return status;

    if ( chanNo > p->fwdfifo->Wmax )
    {
        DEBUG("Channel overflow (passed %td, max %td)", chanNo, p->fwdfifo->Wmax);
        status = LTFATERR_OVERFLOW;

        for (ltfat_int w = p->fwdfifo->Wmax; w < chanNo; w++)
            memset(out[w], 0, len * sizeof * out[w]);

        chanNo = p->fwdfifo->Wmax;
    }

    if ( len > p->bufLenMax )
    {
        DEBUG("Buffer overflow (passed %td, max %td)", len, p->bufLenMax);
        status = LTFATERR_OVERFLOW;

        for (ltfat_int w = 0; w < chanNo; w++)
            memset(out[w] + p->bufLenMax, 0, (len - p->bufLenMax)*sizeof * out[w]);

        len = p->bufLenMax;
    }

    if (!processorCallback)
        processorCallback = &LTFAT_NAME(default_rtdgtreal_processor_callback);

    // Write new data
    samplesWritten = LTFAT_NAME(rtdgtreal_fifo_write)(p->fwdfifo, in, len, chanNo);

    // While there is new data in the input fifo
    while ( LTFAT_NAME(rtdgtreal_fifo_read)(p->fwdfifo, p->buf) > 0 )
    {
        // Transform
        p->fwdtra((void*)p->fwdplan, p->buf, p->fwdfifo->Wmax,
                  p->fftbufIn);

        // Process
        processorCallback(p->userdata, p->fftbufIn, p->fwdplan->M / 2 + 1,
                          p->fwdfifo->Wmax, p->fftbufOut);

        // Reconstruct
        p->backtra((void*)p->backplan, p->fftbufOut, p->backfifo->Wmax, p->buf);

        // Write (and overlap) to out fifo
        LTFAT_NAME(rtidgtreal_fifo_write)(p->backfifo, p->buf);
    }

    // Read sampples for output
    samplesRead = LTFAT_NAME(rtidgtreal_fifo_read)(p->backfifo, len, chanNo, out);

error:
    if (status != LTFATERR_SUCCESS) return status;
    // These should never occur, it would mean internal error
    if ( samplesWritten != len ) return LTFATERR_OVERFLOW;
    else if ( samplesRead != len ) return LTFATERR_UNDERFLOW;
    return status;
}

LTFAT_API int
LTFAT_NAME(rtdgtreal_processor_done)(LTFAT_NAME(rtdgtreal_processor_state)** p)
{
    LTFAT_NAME(rtdgtreal_processor_state)* pp;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);

    pp = *p;
    if (pp->fwdfifo) LTFAT_NAME(rtdgtreal_fifo_done)(&pp->fwdfifo);
    if (pp->backfifo) LTFAT_NAME(rtidgtreal_fifo_done)(&pp->backfifo);
    if (pp->fwdplan) LTFAT_NAME(rtdgtreal_done)(&pp->fwdplan);
    if (pp->backplan) LTFAT_NAME(rtidgtreal_done)(&pp->backplan);
    LTFAT_SAFEFREEALL(pp->buf, pp->fftbufIn, pp->fftbufOut, pp->inTmp, pp->outTmp );

    if (pp->garbageBinSize)
    {
        for (int ii = 0; ii < pp->garbageBinSize; ii++)
            ltfat_safefree(pp->garbageBin[ii]);

        ltfat_safefree(pp->garbageBin);
    }

    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}

LTFAT_API void
LTFAT_NAME(default_rtdgtreal_processor_callback)(void* UNUSED(userdata),
        const LTFAT_COMPLEX* in, int M2, int W, LTFAT_COMPLEX* out)
{
    memcpy(out, in, W * M2 * sizeof * in);
}
