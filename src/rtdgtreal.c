#include "ltfat.h"
#include "ltfat_types.h"

struct LTFAT_NAME(rtdgtreal_plan)
{
    const LTFAT_REAL* g; //!< Window
    const ltfatInt gl; //!< Window length
    const ltfatInt M; //!< Number of FFT channels
    const rtdgt_phasetype ptype; //!< Phase convention
    LTFAT_REAL* fftBuf; //!< Internal buffer
    const ltfatInt fftBufLen; //!< Internal buffer length
    LTFAT_FFTW(plan) pfft; //!< FFTW plan
};

LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_commoninit)(const LTFAT_REAL* g, const ltfatInt gl,
                                 const ltfatInt M, const rtdgt_phasetype ptype,
                                 const ltfat_transformdirection tradir)
{
    LTFAT_NAME(rtdgtreal_plan)* ret = NULL;
    LTFAT_REAL* gshift = NULL;
    LTFAT_REAL* fftBuf = NULL;
    LTFAT_FFTW(plan) pfft = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");

    ltfatInt M2 = M / 2 + 1;
    ltfatInt fftBufLen = gl > 2 * M2 ? gl : 2 * M2;

    CHECKMEM(gshift = ltfat_malloc(gl * sizeof * g));
    CHECKMEM(fftBuf = ltfat_malloc(fftBufLen * sizeof * fftBuf));
    CHECKMEM(ret =  malloc(sizeof * ret));

    LTFAT_NAME_REAL(fftshift)(g, gl, gshift);

    if (LTFAT_FORWARD == tradir)
        pfft = LTFAT_FFTW(plan_dft_r2c_1d)(M, fftBuf, (LTFAT_COMPLEX*)fftBuf,
                                           FFTW_MEASURE);
    else if (LTFAT_INVERSE == tradir)
        pfft = LTFAT_FFTW(plan_dft_c2r_1d)(M, (LTFAT_COMPLEX*)fftBuf, fftBuf,
                                           FFTW_MEASURE);
    else
        CHECKCANTHAPPEN("Unknown transform direction.");

    CHECKINIT(pfft, "FFTW plan creation failed.");

    LTFAT_NAME(rtdgtreal_plan) ret_local =
    {
        .g = gshift, .gl = gl, .M = M, .ptype = ptype,
        .fftBuf = fftBuf, .fftBufLen = fftBufLen, .pfft = pfft
    };

    memcpy(ret, &ret_local, sizeof * ret);
    return ret;
error:
    if (ret) free(ret);
    if (gshift) ltfat_free(gshift);
    if (fftBuf) ltfat_free(fftBuf);
    if (pfft) LTFAT_FFTW(destroy_plan)(pfft);
    return NULL;
}

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                           const ltfatInt M, const rtdgt_phasetype ptype)
{
    return LTFAT_NAME(rtdgtreal_commoninit)(g, gl, M, ptype, LTFAT_FORWARD);
}

LTFAT_EXTERN LTFAT_NAME(rtidgtreal_plan)*
LTFAT_NAME(rtidgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                            const ltfatInt M, const rtdgt_phasetype ptype)
{
    return LTFAT_NAME(rtdgtreal_commoninit)(g, gl, M, ptype, LTFAT_INVERSE);
}


LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_execute)(const LTFAT_NAME(rtdgtreal_plan)* p,
                              const LTFAT_REAL* f, const ltfatInt W,
                              LTFAT_COMPLEX* c)
{
    ltfatInt M = p->M;
    ltfatInt M2 = M / 2 + 1;
    ltfatInt gl = p->gl;
    LTFAT_REAL* fftBuf = p->fftBuf;

    for (ltfatInt w = 0; w < W; w++)
    {
        const LTFAT_REAL* fchan = f + w * gl;
        LTFAT_COMPLEX* cchan = c + w * M2;

        if (p->g)
            for (ltfatInt ii = 0; ii < gl; ii++)
                fftBuf[ii] = fchan[ii] * p->g[ii];

        if (M > gl)
            memset(fftBuf + gl, 0, (M - gl) * sizeof * fftBuf);

        if (gl > M)
            LTFAT_NAME_REAL(fold_array)(fftBuf, gl, M, fftBuf);

        if (p->ptype == LTFAT_RTDGTPHASE_ZERO)
            LTFAT_NAME_REAL(circshift)(fftBuf, M, -(gl / 2), fftBuf );

        LTFAT_FFTW(execute)(p->pfft);

        memcpy(cchan, fftBuf, M2 * sizeof * c);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_execute)(const LTFAT_NAME(rtidgtreal_plan)* p,
                               const LTFAT_COMPLEX* c, const ltfatInt W,
                               LTFAT_REAL* f)
{
    ltfatInt M = p->M;
    ltfatInt M2 = M / 2 + 1;
    ltfatInt gl = p->gl;
    LTFAT_REAL* fftBuf = p->fftBuf;

    for (ltfatInt w = 0; w < W; w++)
    {
        const LTFAT_COMPLEX* cchan = c + w * M2;
        LTFAT_REAL* fchan = f + w * gl;

        memcpy(fftBuf, cchan, M2 * sizeof * cchan);

        LTFAT_FFTW(execute)(p->pfft);

        if (p->ptype == LTFAT_RTDGTPHASE_ZERO)
            LTFAT_NAME_REAL(circshift)(fftBuf, M, gl / 2, fftBuf );

        if (gl > M)
            LTFAT_NAME_REAL(periodize_array)(fftBuf, M , fftBuf, gl);

        if (p->g)
            for (ltfatInt ii = 0; ii < gl; ii++)
                fftBuf[ii] = fftBuf[ii] * p->g[ii];

        memcpy(fchan, fftBuf, gl * sizeof * fchan);
    }
}


LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_done)(LTFAT_NAME(rtdgtreal_plan)* p)
{
    ltfat_free(p->g);
    ltfat_free(p->fftBuf);
    LTFAT_FFTW(destroy_plan)(p->pfft);
    free(p);
}

LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_done)(LTFAT_NAME(rtidgtreal_plan)* p)
{
    LTFAT_NAME(rtdgtreal_done)(p);
}

/* FWD FIFO */

struct LTFAT_NAME(rtdgtreal_fifo)
{
    ltfatInt gl; //!< Window length
    ltfatInt a; //!< Hop size
    LTFAT_REAL* buf; //!< Ring buffer array
    ltfatInt bufLen; //!< Length of the previous
    ltfatInt readIdx; //!< Read pos.
    ltfatInt writeIdx; //!< Write pos.
    ltfatInt Wmax;
};


LTFAT_EXTERN LTFAT_NAME(rtdgtreal_fifo)*
LTFAT_NAME(rtdgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt procDelay,
                                ltfatInt gl, ltfatInt a, ltfatInt Wmax)
{
    LTFAT_REAL* buf = NULL;
    LTFAT_NAME(rtdgtreal_fifo)* ret = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, fifoLen > 0, "fifoLen must be positive");
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");
    CHECK(LTFATERR_BADARG, procDelay >= gl - 1 , "procDelay must be positive");
    CHECK(LTFATERR_BADARG , fifoLen > gl + 1, "fifoLen must be bugger than gl+1");

    CHECKMEM(buf = ltfat_calloc( Wmax * (fifoLen + 1), sizeof * buf));
    CHECKMEM(ret = malloc(sizeof * ret));

    LTFAT_NAME(rtdgtreal_fifo) retloc = {.buf = buf, .bufLen = fifoLen + 1,
                                         .a = a, .gl = gl,
                                         .readIdx =  fifoLen + 1 - (procDelay),
                                         .writeIdx = 0, .Wmax = Wmax
                                        };

    memcpy(ret, &retloc, sizeof * ret);
    return ret;
error:
    if (buf) ltfat_free(buf);
    if (ret) ltfat_free(ret);
    return NULL;
}

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_fifo_done)(LTFAT_NAME(rtdgtreal_fifo)* p)
{
    ltfat_free(p->buf);
    free(p);
}

LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_write)(LTFAT_NAME(rtdgtreal_fifo)* p,
                                 const LTFAT_REAL** buf,
                                 const ltfatInt bufLen, const ltfatInt W)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, bufLen > 0, "bufLen must be positive.");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive.");

    ltfatInt freeSpace = p->readIdx - p->writeIdx - 1;
    if (freeSpace < 0) freeSpace += p->bufLen;

    // Shortcut to the exit
    if (freeSpace == 0) return 0;

    ltfatInt Wact = p->Wmax < W ? p->Wmax : W;

    ltfatInt toWrite = bufLen > freeSpace ? freeSpace : bufLen;
    ltfatInt valid = toWrite;
    ltfatInt over = 0;

    ltfatInt endWriteIdx = p->writeIdx + toWrite;

    if (endWriteIdx > p->bufLen)
    {
        valid = p->bufLen - p->writeIdx;
        over = endWriteIdx - p->bufLen;
    }

    if (valid > 0)
    {
        for (ltfatInt w = 0; w < p->Wmax; w++)
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
        for (ltfatInt w = 0; w < Wact; w++)
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
    return LTFATERR_FAILED;
}

LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_read)(LTFAT_NAME(rtdgtreal_fifo)* p, LTFAT_REAL* buf)
{
    ltfatInt available = p->writeIdx - p->readIdx;
    if (available < 0) available += p->bufLen;

    // Not enough to fill the output buffer
    if (available < p->gl ) return 0;
    ltfatInt toRead = p->gl;

    ltfatInt valid = toRead;
    ltfatInt over = 0;

    ltfatInt endReadIdx = p->readIdx + valid;

    if (endReadIdx > p->bufLen)
    {
        valid = p->bufLen - p->readIdx;
        over = endReadIdx - p->bufLen;
    }

    if (valid > 0)
    {
        for (ltfatInt w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen + p->readIdx;
            memcpy(buf + w * p->gl, pbufchan, valid * sizeof * p->buf );
        }
    }
    if (over > 0)
    {
        for (ltfatInt w = 0; w < p->Wmax; w++)
        {
            memcpy(buf + valid + w * p->gl, p->buf + w * p->bufLen, over * sizeof * p->buf);
        }
    }

    // Only advance by a
    p->readIdx = ( p->readIdx + p->a ) % p->bufLen;

    return toRead;
}

/* BACK FIFO */

struct LTFAT_NAME(rtidgtreal_fifo)
{
    ltfatInt gl; //!< Window length
    ltfatInt a; //!< Hop size
    LTFAT_REAL* buf; //!< Ring buffer array
    ltfatInt bufLen; //!< Length of the previous
    ltfatInt readIdx; //!< Read pos.
    ltfatInt writeIdx; //!< Write pos.
    ltfatInt Wmax;
};


LTFAT_EXTERN LTFAT_NAME(rtidgtreal_fifo)*
LTFAT_NAME(rtidgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl,
                                 ltfatInt a, ltfatInt Wmax)
{
    LTFAT_REAL* buf = NULL;
    LTFAT_NAME(rtidgtreal_fifo)* ret = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, fifoLen > 0, "fifoLen must be positive");
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");
    CHECK(LTFATERR_BADARG , fifoLen > gl + 1, "fifoLen must be bugger than gl+1");


    CHECKMEM(buf = ltfat_calloc( Wmax * (fifoLen + gl + 1), sizeof * buf));
    CHECKMEM(ret = malloc(sizeof * ret));

    LTFAT_NAME(rtidgtreal_fifo) retloc = {.buf = buf, .bufLen = fifoLen + gl + 1,
                                          .a = a, .gl = gl,
                                          .readIdx = 0,
                                          .writeIdx = 0, .Wmax = Wmax
                                         };

    memcpy(ret, &retloc, sizeof * ret);
    return ret;
error:
    if (buf) ltfat_free(buf);
    if (ret) free(ret);
    return NULL;
}

LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_fifo_done)(LTFAT_NAME(rtidgtreal_fifo)* p)
{
    ltfat_free(p->buf);
    free(p);
}

LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_write)(LTFAT_NAME(rtidgtreal_fifo)* p,
                                  const LTFAT_REAL* buf)
{
    ltfatInt freeSpace = p->readIdx - p->writeIdx - 1;
    if (freeSpace < 0) freeSpace += p->bufLen;

    if (freeSpace < p->gl) return 0;

    ltfatInt toWrite = p->gl;
    ltfatInt valid = toWrite;
    ltfatInt over = 0;

    ltfatInt endWriteIdx = p->writeIdx + toWrite;

    if (endWriteIdx > p->bufLen)
    {
        valid = p->bufLen - p->writeIdx;
        over = endWriteIdx - p->bufLen;
    }

    if (valid > 0)
    {
        for (ltfatInt w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + p->writeIdx + w * p->bufLen;
            const LTFAT_REAL* bufchan = buf + w * p->gl;
            for (ltfatInt ii = 0; ii < valid; ii++)
                pbufchan[ii] += bufchan[ii];
        }
    }
    if (over > 0)
    {
        for (ltfatInt w = 0; w < p->Wmax; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen;
            const LTFAT_REAL* bufchan = buf + valid + w * p->gl;
            for (ltfatInt ii = 0; ii < over; ii++)
                pbufchan[ii] += bufchan[ii];
        }
    }

    p->writeIdx = ( p->writeIdx + p->a ) % p->bufLen;

    return toWrite;
}

LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_read)(LTFAT_NAME(rtidgtreal_fifo)* p,
                                 const ltfatInt bufLen, const ltfatInt W,
                                 LTFAT_REAL** buf)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, bufLen > 0, "bufLen must be positive.");
    CHECK(LTFATERR_NOTPOSARG, W > 0, "W must be positive.");

    ltfatInt available = p->writeIdx - p->readIdx;
    if (available < 0) available += p->bufLen;

    ltfatInt toRead = available < bufLen ? available : bufLen;

    ltfatInt valid = toRead;
    ltfatInt over = 0;

    ltfatInt endReadIdx = p->readIdx + valid;

    if (endReadIdx > p->bufLen)
    {
        valid = p->bufLen - p->readIdx;
        over = endReadIdx - p->bufLen;
    }

    // Set the just read samples to zero so that the values are not used in
    // write again
    if (valid > 0)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + p->readIdx + w * p->bufLen;
            memcpy(buf[w], pbufchan, valid * sizeof * p->buf);
            memset(pbufchan, 0, valid * sizeof * p->buf);
        }
    }
    if (over > 0)
    {
        for (ltfatInt w = 0; w < W; w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen;
            memcpy(buf[w] + valid, pbufchan, over * sizeof * p->buf);
            memset(pbufchan, 0, over * sizeof * p->buf);
        }
    }

    p->readIdx = ( p->readIdx + toRead ) % p->bufLen;

    return toRead;
error:
    return LTFATERR_FAILED;
}

/* DGTREAL processor */
struct LTFAT_NAME(rtdgtreal_processor)
{
    LTFAT_NAME(rtdgtreal_processor_callback)
    processorCallback; //!< Custom processor callback
    void* userdata; //!< Callback data
    LTFAT_NAME(rtdgtreal_fifo)* fwdfifo;
    LTFAT_NAME(rtidgtreal_fifo)* backfifo;
    LTFAT_NAME(rtdgtreal_plan)* fwdplan;
    LTFAT_NAME(rtidgtreal_plan)* backplan;
    LTFAT_REAL* buf;
    LTFAT_COMPLEX* fftbufIn;
    LTFAT_COMPLEX* fftbufOut;
    void** garbageBin;
    int garbageBinSize;
};


LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_init)(const LTFAT_REAL* ga, const ltfatInt gal,
                                     const LTFAT_REAL* gs, const ltfatInt gsl,
                                     const ltfatInt a, const ltfatInt M, const ltfatInt Wmax,
                                     LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                     void* userdata)
{
    LTFAT_REAL* buf = NULL;
    LTFAT_COMPLEX* fftbufIn = NULL;
    LTFAT_COMPLEX* fftbufOut = NULL;
    LTFAT_NAME(rtdgtreal_fifo)* fwdfifo = NULL;
    LTFAT_NAME(rtidgtreal_fifo)* backfifo = NULL;
    LTFAT_NAME(rtdgtreal_plan)* fwdplan = NULL;
    LTFAT_NAME(rtidgtreal_plan)* backplan = NULL;
    LTFAT_NAME(rtdgtreal_processor)* ret = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, gal > 0, "gla must be positive");
    CHECK(LTFATERR_NOTPOSARG, gsl > 0, "gls must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");

    CHECKMEM( fftbufIn = ltfat_malloc( Wmax * (M / 2 + 1) * sizeof * fftbufIn));
    CHECKMEM( fftbufOut = ltfat_malloc( Wmax * (M / 2 + 1) * sizeof * fftbufOut));
    CHECKMEM( buf = ltfat_malloc( Wmax * gal * sizeof * buf));
    CHECKINIT( fwdfifo = LTFAT_NAME(rtdgtreal_fifo_init)(11 * gal,
                         gal > gsl ? gal - 1 : gsl - 1 , gal, a, Wmax), "fwd fifo init failed");
    CHECKINIT( backfifo = LTFAT_NAME(rtidgtreal_fifo_init)(11 * gsl, gsl, a, Wmax),
               "back fifo init failed");
    CHECKINIT( fwdplan = LTFAT_NAME(rtdgtreal_init)(ga, gal, M,
                         LTFAT_RTDGTPHASE_ZERO), "fwd plan failed");
    CHECKINIT( backplan = LTFAT_NAME(rtidgtreal_init)(gs, gsl, M,
                          LTFAT_RTDGTPHASE_ZERO), "back plan failed");
    CHECKMEM( ret = malloc(sizeof * ret));

    LTFAT_NAME(rtdgtreal_processor) retLoc =
    {
        .processorCallback = callback, .userdata = userdata,
        .fwdfifo = fwdfifo, .backfifo = backfifo,
        .fwdplan = fwdplan, .backplan = backplan,
        .buf = buf, .fftbufIn = fftbufIn, .fftbufOut = fftbufOut,
        .garbageBin = NULL, .garbageBinSize = 0
    };

    memcpy(ret, &retLoc, sizeof * ret);
    return ret;
error:
    /* if (buf) ltfat_free(buf); */
    /* if (fftbufIn) ltfat_free(fftbufIn); */
    /* if (fftbufOut) ltfat_free(fftbufOut); */
    LTFAT_SAFEFREEALL(buf, fftbufIn, fftbufOut);
    if (ret) free(ret);
    if (fwdfifo) LTFAT_NAME(rtdgtreal_fifo_done)(fwdfifo);
    if (backfifo) LTFAT_NAME(rtidgtreal_fifo_done)(backfifo);
    if (fwdplan) LTFAT_NAME(rtdgtreal_done)(fwdplan);
    if (backplan) LTFAT_NAME(rtidgtreal_done)(backplan);
    return NULL;
}

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_wininit)(LTFAT_FIRWIN win,
                                        const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                        const ltfatInt Wmax,
                                        LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                        void* userdata)
{
    // Enven though the function does not return a status, the check routines require it
    LTFAT_REAL* g = NULL;
    LTFAT_REAL* gd = NULL;
    LTFAT_NAME(rtdgtreal_processor)* ret = NULL;
    void** garbageBin = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0,  "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0,  "M must be positive");
    CHECK(LTFATERR_NOTPOSARG, Wmax > 0, "Wmax must be positive");

    CHECKMEM(g = ltfat_malloc(gl * sizeof * g));
    CHECKMEM(gd = ltfat_malloc(gl * sizeof * gd));
    CHECKMEM(garbageBin = malloc(2 * sizeof(void*)));

    CHECKSTATUS(LTFAT_NAME_REAL(firwin)(win, gl, g), "Call to firwin failed");
    CHECKSTATUS(LTFAT_NAME_REAL(gabdual_painless)(g, gl, a, M, gd),
                "Call to gabdual_painless failed");

    CHECKINIT(ret = LTFAT_NAME(rtdgtreal_processor_init)(g, gl, gd, gl, a, M, Wmax,
                    callback, userdata), "processor_init failed");

    ret->garbageBinSize = 2;
    ret->garbageBin = garbageBin;
    ret->garbageBin[0] = g;
    ret->garbageBin[1] = gd;

    return ret;
error:
    if (g) ltfat_free(g);
    if (gd) ltfat_free(gd);
    if (garbageBin) free(garbageBin);
    // Also status is now set to the proper value
    return NULL;
}

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_processor_setcallback)(LTFAT_NAME(rtdgtreal_processor)* p,
        LTFAT_NAME(rtdgtreal_processor_callback) callback,
        void* userdata)
{
    p->processorCallback = callback;
    p->userdata = userdata;
}

LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_processor_execute)(LTFAT_NAME(rtdgtreal_processor)* p,
                                        const LTFAT_REAL** in,
                                        const ltfatInt len, const ltfatInt chanNo,
                                        LTFAT_REAL** out)
{
    // Get default processor if none was set
    LTFAT_NAME(rtdgtreal_processor_callback) processorCallback =
        p->processorCallback;
    if (!processorCallback)
        processorCallback = &LTFAT_NAME(default_rtdgtreal_processor_callback);

    // Write new data
    ltfatInt samplesWritten = LTFAT_NAME(rtdgtreal_fifo_write)(p->fwdfifo, in, len,
                              chanNo);

    // While there is new data in the input fifo
    while ( LTFAT_NAME(rtdgtreal_fifo_read)(p->fwdfifo, p->buf) )
    {
        // Transform
        LTFAT_NAME(rtdgtreal_execute)(p->fwdplan, p->buf, p->fwdfifo->Wmax,
                                      p->fftbufIn);

        // Process
        (*processorCallback)(p->userdata, p->fftbufIn, p->fwdplan->M / 2 + 1,
                             p->fwdfifo->Wmax, p->fftbufOut);

        // Reconstruct
        LTFAT_NAME(rtidgtreal_execute)(p->backplan, p->fftbufOut, p->backfifo->Wmax,
                                       p->buf);

        // Write (and overlap) to out fifo
        LTFAT_NAME(rtidgtreal_fifo_write)(p->backfifo, p->buf);
    }

    // Read sampples for output
    ltfatInt samplesRead = LTFAT_NAME(rtidgtreal_fifo_read)(p->backfifo, len,
                           chanNo, out);

    if ( samplesWritten != len ) return LTFATERR_OVERFLOW;
    else if ( samplesRead != len ) return LTFATERR_UNDERFLOW;
    return LTFATERR_SUCCESS;
}

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_processor_done)(LTFAT_NAME(rtdgtreal_processor)* p)
{
    LTFAT_NAME(rtdgtreal_fifo_done)(p->fwdfifo);
    LTFAT_NAME(rtidgtreal_fifo_done)(p->backfifo);
    LTFAT_NAME(rtdgtreal_done)(p->fwdplan);
    LTFAT_NAME(rtidgtreal_done)(p->backplan);
    ltfat_free(p->buf);
    ltfat_free(p->fftbufIn);
    ltfat_free(p->fftbufOut);

    if (p->garbageBinSize)
    {
        for (int ii = 0; ii < p->garbageBinSize; ii++)
            ltfat_free(p->garbageBin[ii]);

        free(p->garbageBin);
    }

    free(p);
}

LTFAT_EXTERN void
LTFAT_NAME(default_rtdgtreal_processor_callback)(void* UNUSED(userdata),
        const LTFAT_COMPLEX* in, const int M2, const int W, LTFAT_COMPLEX* out)
{
    memcpy(out, in, W * M2 * sizeof * in);
}
