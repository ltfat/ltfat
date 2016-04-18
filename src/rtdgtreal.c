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
LTFAT_NAME(rtdgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl,
                                ltfatInt a, ltfatInt Wmax)
{
    LTFAT_REAL* buf = ltfat_calloc( Wmax*(fifoLen + 1), sizeof * buf);

    // TODO: This is the worst case delay. 
    // There might be cases for which the following is too pesimistic,
    // but I leave it for now.
    ltfatInt procDelay = gl - 1;
    LTFAT_NAME(rtdgtreal_fifo) retloc = {.buf = buf, .bufLen = fifoLen + 1,
                                         .a = a, .gl = gl,
                                         .readIdx =  fifoLen + 1 - (procDelay),
                                         .writeIdx = 0, .Wmax = Wmax
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
                                 const LTFAT_REAL** buf,
                                 const ltfatInt bufLen, const ltfatInt W)
{
    if (bufLen <= 0) return 0;

    ltfatInt freeSpace = p->readIdx - p->writeIdx - 1;
    if (freeSpace < 0) freeSpace += p->bufLen;

    if (freeSpace == 0) return 0;

    ltfatInt Wact = p->Wmax < W ? p->Wmax: W;

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
        for(ltfatInt w=0;w<p->Wmax;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen + p->writeIdx;
            if(w < Wact)
                memcpy(pbufchan, buf[w], valid * sizeof * p->buf );
            else
                memset(pbufchan, 0, valid * sizeof * p->buf );
        }
    }
    if (over > 0)
    {
        for(ltfatInt w=0;w<Wact;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen;
            if(w < Wact)
                memcpy(pbufchan, buf[w] + valid, over * sizeof * p->buf);
            else
                memset(pbufchan, 0,  over * sizeof * p->buf);
        }
    }
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
    {
        for(ltfatInt w=0;w<p->Wmax;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w * p->bufLen + p->readIdx;
            memcpy(buf + w*p->gl, pbufchan, valid * sizeof * p->buf );
        }
    }
    if (over > 0)
    {
        for(ltfatInt w=0;w<p->Wmax;w++)
        {
           memcpy(buf + valid + w*p->gl, p->buf + w*p->bufLen, over * sizeof * p->buf);
        }
    }

    // Only advance by a
    p->readIdx = ( p->readIdx + p->a ) % p->bufLen;

    return toRead;
}

/* BACK FIFO */

LTFAT_EXTERN LTFAT_NAME(rtidgtreal_fifo)*
LTFAT_NAME(rtidgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl, 
                                 ltfatInt a, ltfatInt Wmax)
{
    // assert(fifoLen > gl)
    // assert(gl > a )

    LTFAT_REAL* buf = ltfat_calloc( Wmax*(fifoLen + gl + 1), sizeof * buf);

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
        for(ltfatInt w=0;w<p->Wmax;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + p->writeIdx + w*p->bufLen;
            const LTFAT_REAL* bufchan = buf + w*p->gl;
            for (ltfatInt ii = 0; ii < valid; ii++)
                pbufchan[ii] += bufchan[ii];
        }
    }
    if (over > 0)
    {
        for(ltfatInt w=0;w<p->Wmax;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w*p->bufLen;
            const LTFAT_REAL* bufchan = buf + valid + w*p->gl;
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
        for(ltfatInt w=0;w<W;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + p->readIdx + w*p->bufLen;
            memcpy(buf[w], pbufchan, valid * sizeof * p->buf);
            memset(pbufchan, 0, valid * sizeof * p->buf);
        }
    }
    if (over > 0)
    {
        for(ltfatInt w=0;w<W;w++)
        {
            LTFAT_REAL* pbufchan = p->buf + w*p->bufLen;
            memcpy(buf[w] + valid, pbufchan, over * sizeof * p->buf);
            memset(pbufchan, 0, over * sizeof * p->buf);
        }
    }

    p->readIdx = ( p->readIdx + toRead ) % p->bufLen;

    return toRead;
}

/* DGTREAL processor */

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_init)(const LTFAT_REAL* g, const LTFAT_REAL* gd,
                                     const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                     const ltfatInt Wmax,
                                     LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                     void* userdata)
{
    LTFAT_NAME(rtdgtreal_processor) retLoc =
    {
        .processorCallback = callback, .userdata = userdata,
        .fwdfifo = LTFAT_NAME(rtdgtreal_fifo_init)(10 * gl, gl, a, Wmax),
        .backfifo = LTFAT_NAME(rtidgtreal_fifo_init)(10 * gl, gl, a, Wmax),
        .fwdplan = LTFAT_NAME(rtdgtreal_init)(g, gl, M, LTFAT_RTDGTPHASE_ZERO),
        .backplan = LTFAT_NAME(rtidgtreal_init)(gd, gl, M, LTFAT_RTDGTPHASE_ZERO),
        .buf = ltfat_malloc( Wmax * gl * sizeof (LTFAT_REAL)),
        .fftbufIn = ltfat_malloc( Wmax * (M / 2 + 1) * sizeof (LTFAT_COMPLEX)),
        .fftbufOut = ltfat_malloc( Wmax * (M / 2 + 1) * sizeof (LTFAT_COMPLEX)),
        .garbageBinSize = 0
    };
    LTFAT_NAME(rtdgtreal_processor)* ret = malloc(sizeof * ret);
    memcpy(ret, &retLoc, sizeof * ret);
    return ret;
}

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_wininit)(LTFAT_FIRWIN win,
                                       const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                       const ltfatInt Wmax,
                                       LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                       void* userdata)
{
    // assert is frame
    // assert is painless
    LTFAT_REAL* g = ltfat_malloc(gl*sizeof*g);
    LTFAT_REAL* gd = ltfat_malloc(gl*sizeof*gd);

    LTFAT_NAME_REAL(firwin)(win,gl,g);
    LTFAT_NAME_REAL(gabdual_painless)(g,gl,a,M,gd);

    LTFAT_NAME(rtdgtreal_processor)* ret =
        LTFAT_NAME(rtdgtreal_processor_init)(g,gd,gl,a,M,Wmax,callback,userdata);

    ret->garbageBinSize = 2;
    ret->garbageBin = malloc(2*sizeof(void*));
    ret->garbageBin[0] = g;
    ret->garbageBin[1] = gd;

    return ret;
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
    ltfatInt samplesWritten = LTFAT_NAME(rtdgtreal_fifo_write)(p->fwdfifo, in, len, chanNo);

    // While there is new data in the input fifo
    while ( LTFAT_NAME(rtdgtreal_fifo_read)(p->fwdfifo, p->buf) )
    {
        // Transform
        LTFAT_NAME(rtdgtreal_execute)(p->fwdplan, p->buf, p->fwdfifo->Wmax, p->fftbufIn);

        // Process
        (*processorCallback)(p->userdata, p->fftbufIn, p->fwdplan->M / 2 + 1,
                             p->fwdfifo->Wmax, p->fftbufOut);

        // Reconstruct
        LTFAT_NAME(rtidgtreal_execute)(p->backplan, p->fftbufOut, p->backfifo->Wmax, p->buf);

        // Write (and overlap) to out fifo
        LTFAT_NAME(rtidgtreal_fifo_write)(p->backfifo, p->buf);
    }

    // Read sampples for output
    ltfatInt samplesRead = LTFAT_NAME(rtidgtreal_fifo_read)(p->backfifo, len, chanNo, out);

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

    if(p->garbageBinSize)
    {
        for(int ii=0;ii<p->garbageBinSize;ii++)
            ltfat_free(p->garbageBin[ii]);

        free(p->garbageBin);
    }

    free(p);
}

LTFAT_EXTERN void
LTFAT_NAME(default_rtdgtreal_processor_callback)(void* UNUSED(userdata),
        const LTFAT_COMPLEX* in, const int M2, LTFAT_COMPLEX* out)
{
    memcpy(out, in, M2 * sizeof * in);
}
