#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_commoninit)(const LTFAT_REAL* g, const ltfatInt gl,
                                 const ltfatInt M, const rtdgt_phasetype ptype,
                                 const ltfat_transformdirection tradir)
{
    LTFAT_NAME(rtdgtreal_plan)* ret = NULL;
    LTFAT_REAL* gshift = NULL;
    LTFAT_REAL* fftBuf = NULL;
    LTFAT_FFTW(plan) pfft = NULL;

    if (!(gshift = ltfat_malloc(gl * sizeof * g))) goto error;

    LTFAT_NAME_REAL(fftshift)(g, gl, gshift);

    ltfatInt M2 = M / 2 + 1;
    ltfatInt fftBufLen = gl > 2 * M2 ? gl : 2 * M2;

    if (!(fftBuf = ltfat_malloc(fftBufLen * sizeof * fftBuf))) goto error;

    if (LTFAT_FORWARD == tradir)
        pfft = LTFAT_FFTW(plan_dft_r2c_1d)(M, fftBuf, (LTFAT_COMPLEX*)fftBuf,
                                           FFTW_MEASURE);
    else if (LTFAT_INVERSE == tradir)
        pfft = LTFAT_FFTW(plan_dft_c2r_1d)(M, (LTFAT_COMPLEX*)fftBuf, fftBuf,
                                           FFTW_MEASURE);
    else
        goto error;

    if (!pfft) goto error;

    LTFAT_NAME(rtdgtreal_plan) ret_local =
    {
        .g = gshift, .gl = gl, .M = M, .ptype = ptype,
        .fftBuf = fftBuf, .fftBufLen = fftBufLen, .pfft = pfft
    };

    if (!(ret = malloc(sizeof * ret))) goto error;
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

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_fifo)*
LTFAT_NAME(rtdgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl, ltfatInt a)
{
    LTFAT_REAL* buf = ltfat_calloc( (fifoLen + 1), sizeof * buf);
 
    // TODO: This is the worst case delay. 
    // There might be cases for which the following is too pesimistic,
    // but I leave it for now.
    ltfatInt procDelay = gl - 1;
    LTFAT_NAME(rtdgtreal_fifo) retloc = {.buf = buf, .bufLen = fifoLen + 1,
                                         .a = a, .gl = gl,
                                         .readIdx =  fifoLen + 1 - (procDelay),
                                         .writeIdx = 0
                                        };

    LTFAT_NAME(rtdgtreal_fifo)* ret = malloc(sizeof * ret);
    memcpy(ret, &retloc, sizeof * ret);
    return ret;
}

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_fifo_done)(LTFAT_NAME(rtdgtreal_fifo)* p)
{
    ltfat_free(p->buf);
    free(p);
}

LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_write)(LTFAT_NAME(rtdgtreal_fifo)* p,
                                 const ltfatInt bufLen, const LTFAT_REAL* buf)
{
    if (bufLen <= 0) return 0;

    ltfatInt freeSpace = p->readIdx - p->writeIdx - 1;
    if (freeSpace < 0) freeSpace += p->bufLen;

    if (freeSpace == 0) return 0;

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
        memcpy(p->buf + p->writeIdx, buf, valid * sizeof * p->buf );
    if (over > 0)
        memcpy(p->buf, buf + valid, over * sizeof * p->buf);

    p->writeIdx = ( p->writeIdx + toWrite ) % p->bufLen;

    return toWrite;
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
        memcpy(buf, p->buf + p->readIdx, valid * sizeof * p->buf );
    if (over > 0)
        memcpy(buf + valid, p->buf, over * sizeof * p->buf);

    // Only advance by a
    p->readIdx = ( p->readIdx + p->a ) % p->bufLen;

    return toRead;
}

/* BACK FIFO */

LTFAT_EXTERN LTFAT_NAME(rtidgtreal_fifo)*
LTFAT_NAME(rtidgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl, ltfatInt a)
{
    // assert(fifoLen > gl)
    // assert(gl > a )

    LTFAT_REAL* buf = ltfat_calloc((fifoLen + gl + 1), sizeof * buf);

    // If gl is not integer divisible by a, we have to work with the worst case
    // delay
    /* ltfatInt procDelay = gl % a ? gl - 1 : gl - a; */

    // Initialize readIdx and writeIdx one slot apart..
    LTFAT_NAME(rtidgtreal_fifo) retloc = {.buf = buf, .bufLen = fifoLen + gl + 1,
                                          .a = a, .gl = gl,
                                          .readIdx = 0,
                                          .writeIdx = 0
                                         };

    LTFAT_NAME(rtidgtreal_fifo)* ret = malloc(sizeof * ret);
    memcpy(ret, &retloc, sizeof * ret);
    return ret;
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
        LTFAT_REAL* pbuf = p->buf + p->writeIdx;
        for (ltfatInt ii = 0; ii < valid; ii++)
            pbuf[ii] += buf[ii];
    }
    if (over > 0)
    {
        const LTFAT_REAL* vbuf = buf + valid;
        for (ltfatInt ii = 0; ii < over; ii++)
            p->buf[ii] += vbuf[ii];
    }

    p->writeIdx = ( p->writeIdx + p->a ) % p->bufLen;

    return toWrite;
}

LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_read)(LTFAT_NAME(rtidgtreal_fifo)* p,
                                 const ltfatInt bufLen, LTFAT_REAL* buf)
{
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
        memcpy(buf, p->buf + p->readIdx, valid * sizeof * p->buf);
        memset(p->buf + p->readIdx, 0, valid * sizeof * p->buf);
    }
    if (over > 0)
    {
        memcpy(buf + valid, p->buf, over * sizeof * p->buf);
        memset(p->buf, 0, over * sizeof * p->buf);
    }

    p->readIdx = ( p->readIdx + toRead ) % p->bufLen;

    return toRead;
}

/* DGTREAL processor */

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_init)(const LTFAT_REAL* g, const LTFAT_REAL* gd,
                                     const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                     LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                     void* userdata)
{
    LTFAT_NAME(rtdgtreal_processor) retLoc =
    {
        .processorCallback = callback, .userdata = userdata,
        .fwdfifo = LTFAT_NAME(rtdgtreal_fifo_init)(10 * gl, gl, a),
        .backfifo = LTFAT_NAME(rtidgtreal_fifo_init)(10 * gl, gl, a),
        .fwdplan = LTFAT_NAME(rtdgtreal_init)(g, gl, M, LTFAT_RTDGTPHASE_ZERO),
        .backplan = LTFAT_NAME(rtidgtreal_init)(gd, gl, M, LTFAT_RTDGTPHASE_ZERO),
        .buf = ltfat_malloc(gl * sizeof (LTFAT_REAL)),
        .fftbufIn = ltfat_malloc((M / 2 + 1) * sizeof (LTFAT_COMPLEX)),
        .fftbufOut = ltfat_malloc((M / 2 + 1) * sizeof (LTFAT_COMPLEX))
    };
    LTFAT_NAME(rtdgtreal_processor)* ret = malloc(sizeof * ret);
    memcpy(ret, &retLoc, sizeof * ret);
    return ret;
}


LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_processor_execute)(LTFAT_NAME(rtdgtreal_processor)* p,
                                        const LTFAT_REAL* in, const ltfatInt len, LTFAT_REAL* out)
{
    // Get default processor if none was set
    LTFAT_NAME(rtdgtreal_processor_callback) processorCallback =
        p->processorCallback;
    if (!processorCallback)
        processorCallback = &LTFAT_NAME(default_rtdgtreal_processor_callback);

    // Write new data
    ltfatInt samplesWritten = LTFAT_NAME(rtdgtreal_fifo_write)(p->fwdfifo, len, in);

    // While there is new data in the input fifo
    while ( LTFAT_NAME(rtdgtreal_fifo_read)(p->fwdfifo, p->buf) )
    {
        // Transform
        LTFAT_NAME(rtdgtreal_execute)(p->fwdplan, p->buf, 1, p->fftbufIn);

        // Process
        (*processorCallback)(p->userdata, p->fftbufIn, p->fwdplan->M / 2 + 1,
                             p->fftbufOut);

        // Reconstruct
        LTFAT_NAME(rtidgtreal_execute)(p->backplan, p->fftbufOut, 1, p->buf);

        // Write (and overlap) to out fifo
        LTFAT_NAME(rtidgtreal_fifo_write)(p->backfifo, p->buf);
    }

    // Read sampples for output
    ltfatInt samplesRead = LTFAT_NAME(rtidgtreal_fifo_read)(p->backfifo, len, out);

    if ( samplesWritten != len ) return LTFAT_RTDGT_STREAMOVERFLOW;
    else if ( samplesRead != len ) return LTFAT_RTDGT_STREAMUNDERFLOW;
    return LTFAT_RTDGT_STREAMOK;
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
    free(p);
}

LTFAT_EXTERN void
LTFAT_NAME(default_rtdgtreal_processor_callback)(void* UNUSED(userdata),
        const LTFAT_COMPLEX* in, const ltfatInt M2, LTFAT_COMPLEX* out)
{
    memcpy(out, in, M2 * sizeof * in);
}
