#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealmp_private.h"
#include "slicingbuf_private.h"
#include "slidgtrealmp_private.h"

LTFAT_API int
LTFAT_NAME(slidgtrealmp_init)(
    LTFAT_NAME(dgtrealmp_state)* mpstate,
    LTFAT_NAME(slicing_processor_state)* slistate,
    LTFAT_NAME(slidgtrealmp_state)** pout)
{
    int status = LTFATERR_FAILED;

    LTFAT_NAME(slidgtrealmp_state)* p = NULL;

    CHECKNULL(mpstate); CHECKNULL(slistate); CHECKNULL(pout);
    CHECKMEM( p = LTFAT_NEW(LTFAT_NAME(slidgtrealmp_state)));
    CHECKMEM( p->couttmp = LTFAT_NEWARRAY(LTFAT_COMPLEX*, mpstate->P));

    for(ltfat_int pidx = 0; pidx< mpstate->P; pidx++)
        CHECKMEM( p->couttmp[pidx] = LTFAT_NAME_COMPLEX(malloc)(
                    mpstate->M2[pidx]*mpstate->L/mpstate->a[pidx]));

    LTFAT_NAME(slicing_processor_setcallback)( slistate,
               &LTFAT_NAME(slidgtrealmp_execute_callback), p);

    p->mpstate = mpstate;
    p->slistate = slistate;
    *pout = p;
    return LTFATERR_SUCCESS;
error:
    if(p) LTFAT_NAME(slidgtrealmp_done)(&p);
    return status;
}

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute)(
    LTFAT_NAME(slidgtrealmp_state)* p,
    const LTFAT_REAL* in[], ltfat_int inLen, ltfat_int chanNo,
    LTFAT_REAL* out[])
{
    int status = LTFATERR_FAILED;
    CHECKNULL(p);
    return LTFAT_NAME(slicing_processor_execute)( p->slistate, in, inLen, chanNo, out);
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute_compact)(
    LTFAT_NAME(slidgtrealmp_state)* p,
    const LTFAT_REAL in[], ltfat_int inLen, ltfat_int chanNo,
    LTFAT_REAL out[])
{
    int status = LTFATERR_FAILED;
    CHECKNULL(p);
    return LTFAT_NAME(slicing_processor_execute_compact)( p->slistate, in, inLen, chanNo, out);
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(slidgtrealmp_done)(LTFAT_NAME(slidgtrealmp_state)** p)
{
    LTFAT_NAME(slidgtrealmp_state)* pp;
    int status = LTFATERR_FAILED;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;
    ltfat_free(pp);
    pp = NULL;
    return LTFATERR_SUCCESS;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(slidgtrealmp_reset)(
        LTFAT_NAME(slidgtrealmp_state)* p);

int
LTFAT_NAME(slidgtrealmp_execute_callback)(void* userdata,
        const LTFAT_REAL in[], int winLen, int UNUSED(taperLen), int UNUSED(zpadLen), int W, LTFAT_REAL out[])
{

    LTFAT_NAME(slidgtrealmp_state)* p =
        (LTFAT_NAME(slidgtrealmp_state)*) userdata;

    for(ltfat_int w = 0; w < W; w++)
    {
        LTFAT_NAME(dgtrealmp_reset)( p->mpstate, in + w*winLen);
        LTFAT_NAME(dgtrealmp_execute)(
                    p->mpstate, in + w*winLen, p->couttmp, out + w*winLen);
    }
    return  0;
}

/* LTFAT_API int */
/* LTFAT_NAME(slidgtrealmp_setnitercallback)( */
/*         LTFAT_NAME(slidgtrealmp_state)* p, */
/*         LTFAT_NAME(slidgtrealmp_niter_callback)* callback, */
/*         void* userdata) */
/* { */
/*     int status = LTFATERR_FAILED; */
/*     CHECKNULL(p); */
/*     p->callback = callback; */
/*     p->userdata = userdata; */
/*     return LTFATERR_SUCCESS; */
/* error: */
/*     return status; */
/* } */




