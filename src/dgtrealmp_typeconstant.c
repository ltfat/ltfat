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

LTFAT_API int
ltfat_dgtrealmp_setpar_hint(ltfat_dgtrealmp_params* params,
                            ltfat_dgtrealmp_hint hint)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    CHECK(LTFATERR_BADARG, ltfat_dgtrealmp_hint_isvalid(hint),
          "Invalid hint passed (passed %d)", hint);

    params->hint = hint;

error:
    return status;
}

LTFAT_API int
ltfat_dgtrealmp_setpar_alg(ltfat_dgtrealmp_params* params,
                           ltfat_dgtrealmp_alg alg)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    CHECK(LTFATERR_BADARG, ltfat_dgtrealmp_alg_isvalid(alg),
          "Invalid hint passed (passed %d)", alg);

    params->alg = alg;
error:
    return status;
}

LTFAT_API int
ltfat_dgtrealmp_setpar_maxatoms(ltfat_dgtrealmp_params* params,
                                size_t maxatoms)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);

    CHECK(LTFATERR_NOTPOSARG, maxatoms > 0, "maxatoms must be greater than 0");
    params->maxatoms = maxatoms;

    if (params->maxit == 0)
        params->maxit = 2 * maxatoms;
error:
    return status;
}

LTFAT_API int
ltfat_dgtrealmp_setpar_errtoldb(ltfat_dgtrealmp_params* params,
                                double errtoldb)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    CHECK(LTFATERR_BADARG, errtoldb <= 0, "errtoldb must be lower than 0");
    params->errtoldb = errtoldb;

error:
    return status;
}

LTFAT_API int
ltfat_dgtrealmp_setpar_iterstep(ltfat_dgtrealmp_params* params,
                                size_t iterstep)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    CHECK(LTFATERR_NOTPOSARG, iterstep > 0, "iterstep must be greater than 0");
    params->iterstep = iterstep;
error:
    return status;
}

LTFAT_API int
ltfat_dgtrealmp_setpar_kernrelthr(ltfat_dgtrealmp_params* params,
                                  double thr)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(params);
    CHECK(LTFATERR_BADARG, thr <= 1 && thr >= 0,
          "Relative threshold must be in range [0-1] (passed %s)", thr);
    params->kernrelthr = thr;
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
    params->errtoldb = -40.0;
    params->kernrelthr = 1e-4;
    params->verbose = 1;
    params->maxatoms = 0;
    params->maxit = 0;
    params->iterstep = 100;
    params->treelevels = 10;
    params->ptype = LTFAT_FREQINV;
error:
    return status;
}

int
ltfat_dgtrealmp_hint_isvalid(ltfat_dgtrealmp_hint in)
{
    int isvalid = 0;

    switch (in)
    {
    case ltfat_dgtrealmp_singlemod:
    case ltfat_dgtrealmp_allmods:
        isvalid = 1;
    }

    return isvalid;
}

int
ltfat_dgtrealmp_alg_isvalid(ltfat_dgtrealmp_alg in)
{
    int isvalid = 0;

    switch (in)
    {
    case ltfat_dgtrealmp_alg_mp:
    case ltfat_dgtrealmp_alg_locomp:
    case ltfat_dgtrealmp_alg_cyclicmp:
        isvalid = 1;
    }

    return isvalid;
}
