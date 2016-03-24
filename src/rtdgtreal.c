#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                           const ltfatInt M, const rtdgt_phasetype ptype)
{
    LTFAT_NAME(rtdgtreal_plan)* ret = malloc(sizeof * ret);

    LTFAT_REAL* gshift = ltfat_malloc(gl * sizeof * g);
    LTFAT_NAME_REAL(fftshift)(g, gl, gshift);

    ltfatInt M2 = M / 2 + 1;
    ltfatInt fftBufLen = gl > 2 * M2 ? gl : 2 * M2;
    LTFAT_REAL* fftBuf = ltfat_malloc(fftBufLen * sizeof * fftBuf);

    LTFAT_FFTW(plan) pfft = LTFAT_FFTW(plan_dft_r2c_1d)(M, fftBuf, fftBuf,
                            FFTW_MEASURE);

    LTFAT_NAME(rtdgtreal_plan)* ret_local =
    {
        .g = gshift, .gl = gl, .M = M, .ptype = ptype,
        .fftBuf = fftBuf, .fftBufLen = bufLen, .pfft = pfft
    };

    memcpy(ret, ret_local, sizeof * ret);
    return ret;
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
        LTFAT_COMPLEX cchan = c + w * M2;

        for (ltfatInt ii = 0; ii < gl; ii++)
            fftBuf[ii] = fchan[ii] * g[ii];

        if (M > gl)
            memset(fftBuf + gl, 0, (M - gl)*sizeof * fftBuf);

        if (gl > M)
            LTFAT_NAME_REAL(fold_array)(fftBuf, gl, M, fftBuf);

        if (LTFAT_RTDGTPHASE_ZERO)
            LTFAT_NAME_REAL(circshift)(fftBuf, M, -(gl / 2), fftBuf );

        LTFAT_FFTW(execute)(p->pfft);

        memcpy(cchan, fftBuf, M2 * sizeof * c);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(irtdgtreal_execute)(const LTFAT_NAME(rtdgtreal_plan)* p,
                               const LTFAT_COMPLEX* c, const ltfatInt W,
                               LTFAT_REAL* f)
{
    ltfatInt M = p->M;
    ltfatInt M2 = M / 2 + 1;
    ltfatInt gl = p->gl;
    LTFAT_COMPLEX* fftBufCmplx = (LTFAT_COMPLEX*) p->fftBuf;
    LTFAT_REAL* fftBuf = p->fftBuf;

    for (ltfatInt w = 0; w < W; w++)
    {
        const LTFAT_COMPLEX cchan = c + w * M2;
        LTFAT_REAL* fchan = f + w * gl;

        memcpy(fftBuf, cchan, M2 * sizeof(cchan));

        LTFAT_FFTW(execute)(p->pfft);

        if (LTFAT_RTDGTPHASE_ZERO)
            LTFAT_NAME_REAL(circshift)(fftBuf, M, gl / 2, fftBuf );

        if (gl > M)
            LTFAT_NAME_REAL(periodize_array)(fftBuf, M , gl, fftBuf);

        for (ltfatInt ii = 0; ii < gl; ii++)
            fftBuf[ii] = fftBuf[ii] * g[ii] / M;

        memcpy(fchan, fftBuf, gl * sizeof * f);
    }


}


LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_done)(const LTFAT_NAME(rtdgtreal_plan)* p)
{
    ltfat_free(p->g);
    ltfat_free(fftBuf);
    LTFAT_FFTW(destroy_plan)(p->pfft);
    free(p);
}
