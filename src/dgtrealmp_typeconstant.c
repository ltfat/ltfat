#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealmp_private.h"

LTFAT_API ltfat_dgtrealmp_params*
ltfat_dgtrealmp_params_allocdef()
{
    ltfat_dgtrealmp_params* params;
    int status = LTFATERR_SUCCESS;
    CHECKMEM( params = LTFAT_NEW(ltfat_dgtrealmp_params));

    ltfat_dgtrealmp_params_defaults(params);

error:
    return params;
}

LTFAT_API int
ltfat_dgtrealmp_params_free(ltfat_dgtrealmp_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    ltfat_free(params);
error:
    return status;
}

int
ltfat_dgtrealmp_params_defaults(ltfat_dgtrealmp_params* params)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    params->hint = ltfat_dgtrealmp_allmods;
    params->alg = ltfat_dgtrealmp_alg_mp;
    params->errtol = -40.0;
    params->kernrelthr = 1e-4;
    params->verbose = 1;
    params->maxatoms = 0;
    params->maxit = 0;
    params->iterstep = 10;
error:
    return status;
}
