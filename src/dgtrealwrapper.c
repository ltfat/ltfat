#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealwrapper_private.h"

int
LTFAT_NAME(idgtreal_long_execute_wrapper)(void* plan,
        const LTFAT_COMPLEX* c, ltfat_int UNUSED(L), ltfat_int UNUSED(W), LTFAT_REAL* f)
{
    return LTFAT_NAME(idgtreal_long_execute_newarray)(
               (LTFAT_NAME(idgtreal_long_plan)*) plan, c, f);
}

int
LTFAT_NAME(dgtreal_long_execute_wrapper)(void* plan,
        const LTFAT_REAL* f, ltfat_int UNUSED(L), ltfat_int UNUSED(W), LTFAT_COMPLEX* c)
{
    return LTFAT_NAME(dgtreal_long_execute_newarray)(
               (LTFAT_NAME(dgtreal_long_plan)*) plan, f, c);
}

int
LTFAT_NAME(idgtreal_fb_execute_wrapper)(void* plan,
                                        const LTFAT_COMPLEX* c, ltfat_int L, ltfat_int W, LTFAT_REAL* f)
{
    return LTFAT_NAME(idgtreal_fb_execute)(
               (LTFAT_NAME(idgtreal_fb_plan)*) plan, c, L, W, f);
}

int
LTFAT_NAME(dgtreal_fb_execute_wrapper)(void* plan,
                                       const LTFAT_REAL* f, ltfat_int L, ltfat_int W, LTFAT_COMPLEX* c)
{
    return LTFAT_NAME(dgtreal_fb_execute)(
               (LTFAT_NAME(dgtreal_fb_plan)*) plan, f, L, W, c);
}

int
LTFAT_NAME(idgtreal_long_done_wrapper)(void** plan)
{
    return LTFAT_NAME(idgtreal_long_done)( (LTFAT_NAME(idgtreal_long_plan)**) plan);
}

int
LTFAT_NAME(dgtreal_long_done_wrapper)(void** plan)
{
    return LTFAT_NAME(dgtreal_long_done)( (LTFAT_NAME(dgtreal_long_plan)**) plan);
}

int
LTFAT_NAME(idgtreal_fb_done_wrapper)(void** plan)
{
    return LTFAT_NAME(idgtreal_fb_done)((LTFAT_NAME(idgtreal_fb_plan)**) plan);
}

int
LTFAT_NAME(dgtreal_fb_done_wrapper)(void** plan)
{
    return LTFAT_NAME(dgtreal_fb_done)((LTFAT_NAME(dgtreal_fb_plan)**) plan);
}

LTFAT_API int
LTFAT_NAME(dgtreal_execute_proj)(
    LTFAT_NAME(dgtreal_plan)* p, const LTFAT_COMPLEX cin[],
    LTFAT_COMPLEX cout[])
{
    int status = LTFATERR_SUCCESS;
    CHECKSTATUS( p->backtra(p->backtra_userdata, cin, p->L, p->W, p->f),
                 "Back transform failed");
    CHECKSTATUS( p->fwdtra(p->fwdtra_userdata, p->f, p->L, p->W,  cout),
                 "Forward transform failed");
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(dgtreal_execute_syn)(
    LTFAT_NAME(dgtreal_plan)* p, const LTFAT_COMPLEX c[], LTFAT_REAL f[])
{
    return p->backtra(p->backtra_userdata, c, p->L, p->W, f);
}

LTFAT_API int
LTFAT_NAME(dgtreal_execute_ana)(
    LTFAT_NAME(dgtreal_plan)* p, const LTFAT_REAL f[], LTFAT_COMPLEX c[])
{
    return p->fwdtra(p->fwdtra_userdata, f, p->L, p->W,  c);
}

LTFAT_API int
LTFAT_NAME(dgtreal_done)(LTFAT_NAME(dgtreal_plan)** p)
{
    int status = LTFATERR_SUCCESS;
    LTFAT_NAME(dgtreal_plan)* pp;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    if (pp->fwdtra_userdata)
        CHECKSTATUS( pp->fwddonefunc(&pp->fwdtra_userdata),
                     "Forward transform done function failed");

    if (pp->backtra_userdata)
        CHECKSTATUS( pp->backdonefunc(&pp->backtra_userdata),
                     "Back trnasform done function failed");

    ltfat_safefree(pp->f);
    ltfat_free(pp);
    pp = NULL;
error:
    return status;
}


LTFAT_API int
LTFAT_NAME(dgtreal_init)(const LTFAT_REAL g[], ltfat_int gl, ltfat_int L,
                         ltfat_int W, ltfat_int a, ltfat_int M,
                         LTFAT_REAL f[], LTFAT_COMPLEX c[],
                         ltfat_dgtreal_params* params, LTFAT_NAME(dgtreal_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    int ispainless = gl <= M;
    LTFAT_REAL* g2 = NULL;
    ltfat_int g2l = 0;

    ltfat_int minL = ltfat_lcm(a, M);

    CHECK(LTFATERR_BADTRALEN, !(L % minL),
          "L must divisible by lcm(a,M)=%d.", minL);

    if (ispainless)
    {
        // The length of the dual window is guaranteed to be gl
        g2l = gl;
        CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(gl));
        CHECKSTATUS( LTFAT_NAME(gabdual_painless)(g, gl, a, M, g2),
                     "Gabdual painless call failed");
    }
    else
    {
#ifndef NOBLASLAPACK
        g2l = L;
        CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(L));
        LTFAT_NAME(fir2long)(g, gl, L, g2);
        CHECKSTATUS( LTFAT_NAME(gabdual_long)(g, L, a, M, g2),
                     "Gabdual long failed");
#else
        CHECK( LTFATERR_NOTSUPPORTED, 0, "Non-painless support was not compiled.");
#endif
    }

    CHECKSTATUS(
        LTFAT_NAME(dgtreal_init_gen)(g, gl, g2, g2l, L, W, a, M, f, c, params, pout),
        "dgtreal_init_gen failed");

    ltfat_free(g2);
    return status;
error:
    ltfat_safefree(g2);
    return status;
}

LTFAT_API int
LTFAT_NAME(dgtreal_init_gen)(const LTFAT_REAL ga[], ltfat_int gal,
                             const LTFAT_REAL gs[], ltfat_int gsl,
                             ltfat_int L, ltfat_int W, ltfat_int a, ltfat_int M,
                             LTFAT_REAL f[], LTFAT_COMPLEX c[],
                             ltfat_dgtreal_params* params, LTFAT_NAME(dgtreal_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    LTFAT_NAME(dgtreal_plan)* p = NULL;
    ltfat_dgtreal_params paramsLoc;

    LTFAT_REAL* g2 = NULL;

    ltfat_int minL = ltfat_lcm(a, M);

    if (params)
        paramsLoc = *params;
    else
        ltfat_dgtreal_params_defaults(&paramsLoc);

    CHECK(LTFATERR_BADTRALEN, !(L % minL),
          "L must divisible by lcm(a,M)=%d.", minL);

    CHECKMEM( p = LTFAT_NEW(LTFAT_NAME(dgtreal_plan)) );
    p->M = M, p->a = a, p->L = L, p->W = W, p->c = c; p->f = f;

    if (ltfat_dgtreal_long == paramsLoc.hint)
    {
        // Make the dual window longer if it is not already
        CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(L) );
        LTFAT_NAME(fir2long)(gs, gsl, L, g2);

        p->backtra = &LTFAT_NAME(idgtreal_long_execute_wrapper);
        p->backdonefunc = &LTFAT_NAME(idgtreal_long_done_wrapper);

        CHECKSTATUS(
            LTFAT_NAME(idgtreal_long_init)(c, g2, L, W, a, M, p->f, paramsLoc.ptype,
                                           paramsLoc.fftw_flags,
                                           (LTFAT_NAME(idgtreal_long_plan)**)&p->backtra_userdata),
            "idgtreal long init failed"
        );

        p->fwdtra = &LTFAT_NAME(dgtreal_long_execute_wrapper);
        p->fwddonefunc = &LTFAT_NAME(dgtreal_long_done_wrapper);

        // Ensure the original window is long enough
        LTFAT_NAME(fir2long)(ga, gal, L, g2);

        CHECKSTATUS(
            LTFAT_NAME(dgtreal_long_init)(p->f, g2, L, W, a, M, c, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(dgtreal_long_plan)**)&p->fwdtra_userdata),
            "dgtreal long init failed");

        ltfat_safefree(g2);
    }
    else if ( ltfat_dgtreal_fb == paramsLoc.hint )
    {
        // Use _fb functions only
        p->backtra = &LTFAT_NAME(idgtreal_fb_execute_wrapper);
        p->backdonefunc = &LTFAT_NAME(idgtreal_fb_done_wrapper);

        CHECKSTATUS(
            LTFAT_NAME(idgtreal_fb_init)( gs, gsl, a, M, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(idgtreal_fb_plan)**)&p->backtra_userdata),
            "idgtreal fb init failed");

        p->fwdtra = &LTFAT_NAME(dgtreal_fb_execute_wrapper);
        p->fwddonefunc = &LTFAT_NAME(dgtreal_fb_done_wrapper);

        CHECKSTATUS(
            LTFAT_NAME(dgtreal_fb_init)( ga, gal, a, M, paramsLoc.ptype,
                                         paramsLoc.fftw_flags,
                                         (LTFAT_NAME(dgtreal_fb_plan)**)&p->fwdtra_userdata),
            "dgtreal fb init failed");

    }
    else if ( ltfat_dgtreal_auto == paramsLoc.hint )
    {
        // Decide whether to use _fb or _long depending on the window lengths
        if (gsl < L)
        {
            p->backtra = &LTFAT_NAME(idgtreal_fb_execute_wrapper);
            p->backdonefunc = &LTFAT_NAME(idgtreal_fb_done_wrapper);

            LTFAT_NAME(idgtreal_fb_init)( gs, gsl, a, M, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(idgtreal_fb_plan)**)&p->backtra_userdata);

        }
        else
        {
            p->backtra = &LTFAT_NAME(idgtreal_long_execute_wrapper);
            p->backdonefunc = &LTFAT_NAME(idgtreal_long_done_wrapper);

            CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(L) );
            LTFAT_NAME(fir2long)(gs, gsl, L, g2);

            LTFAT_NAME(idgtreal_long_init)(c, g2, L, W, a, M, p->f, paramsLoc.ptype,
                                           paramsLoc.fftw_flags,
                                           (LTFAT_NAME(idgtreal_long_plan)**)&p->backtra_userdata);

            ltfat_safefree(g2);
        }

        if (gal < L)
        {
            p->fwdtra = &LTFAT_NAME(dgtreal_fb_execute_wrapper);
            p->fwddonefunc = &LTFAT_NAME(dgtreal_fb_done_wrapper);

            LTFAT_NAME(dgtreal_fb_init)(ga, gal, a, M, paramsLoc.ptype,
                                        paramsLoc.fftw_flags,
                                        (LTFAT_NAME(dgtreal_fb_plan)**)&p->fwdtra_userdata);

        }
        else
        {
            p->fwdtra = &LTFAT_NAME(dgtreal_long_execute_wrapper);
            p->fwddonefunc = &LTFAT_NAME(dgtreal_long_done_wrapper);

            CHECKMEM( g2 = LTFAT_NAME_REAL(malloc)(L) );
            LTFAT_NAME(fir2long)(ga, gal, L, g2);


            LTFAT_NAME(dgtreal_long_init)(p->f, g2, L, W, a, M, c, paramsLoc.ptype,
                                          paramsLoc.fftw_flags,
                                          (LTFAT_NAME(dgtreal_long_plan)**)&p->fwdtra_userdata);

            ltfat_safefree(g2);
        }
    }
    else
    {
        CHECKCANTHAPPEN("No such dgtreal hint");
    }

    *pout = p;

    return status;
error:
    ltfat_safefree(g2);
    if (p) LTFAT_NAME(dgtreal_done)(&p);
    return status;


}
