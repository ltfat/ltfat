#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealwrapper_private.h"

#include "ltfat/thirdparty/fftw3.h"

int
ltfat_dgtreal_params_defaults(ltfat_dgtreal_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->ptype = LTFAT_FREQINV;
    params->fftw_flags = FFTW_ESTIMATE;
    params->hint = ltfat_dgtreal_auto;
error:
    return status;
}

LTFAT_API ltfat_dgtreal_params*
ltfat_dgtreal_params_allocdef()
{
    ltfat_dgtreal_params* params;
    int status = LTFATERR_SUCCESS;
    CHECKMEM( params = LTFAT_NEW(ltfat_dgtreal_params));

    ltfat_dgtreal_params_defaults(params);
error:
    return params;
}

LTFAT_API int
ltfat_dgtreal_setpar_phaseconv(ltfat_dgtreal_params* params,
                                   ltfat_phaseconvention ptype)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->ptype = ptype;
error:
    return status;
}

LTFAT_API int
ltfat_dgtreal_setpar_fftwflags(ltfat_dgtreal_params* params,
                                   unsigned fftw_flags)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->fftw_flags = fftw_flags;
error:
    return status;

}

LTFAT_API int
ltfat_dgtreal_setpar_hint(ltfat_dgtreal_params* params,
                              ltfat_dgtreal_hint hint)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->hint = hint;
error:
    return status;
}

/* LTFAT_API int */
/* ltfat_dgtreal_setpar_normalizewin(ltfat_dgtreal_params* params, */
/*                                       int do_normalize_win) */
/* { */
/*     int status = LTFATERR_SUCCESS; */
/*     CHECKNULL(params); */
/*     params->normalize_win = do_normalize_win; */
/* error: */
/*     return status; */
/* } */

LTFAT_API int
ltfat_dgtreal_params_free(ltfat_dgtreal_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    ltfat_free(params);
error:
    return status;
}
