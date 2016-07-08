#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

LTFAT_EXTERN int
LTFAT_NAME(gabtight_long)(const LTFAT_TYPE* g,
                          const ltfatInt L, const ltfatInt R, const ltfatInt a,
                          const ltfatInt M, LTFAT_TYPE* gd)
{
    LTFAT_COMPLEX* gf = NULL;
    LTFAT_COMPLEX* gdf = NULL;

    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L (passed %d) must be positive.", L);
    CHECK(LTFATERR_NOTPOSARG, R > 0, "R (passed %d) must be positive.", R);
    // a,M, g and gd are checked further

    CHECKMEM( gf = ltfat_malloc(L * R * sizeof * gf));
    CHECKMEM( gdf = ltfat_malloc(L * R * sizeof * gdf));

#ifdef LTFAT_COMPLEXTYPE

    CHECKSTATUS( LTFAT_NAME(wfac)(g, L, R, a, M, gf), "wfac failed");
    LTFAT_NAME_REAL(gabtight_fac)(gf, L, R, a, M, gdf);
    CHECKSTATUS( LTFAT_NAME(iwfac)(gdf, L, R, a, M, gd), "iwfac failed");

#else

    LTFAT_NAME_REAL(wfacreal)(g, L, R, a, M, gf);
    LTFAT_NAME_REAL(gabtightreal_fac)((const LTFAT_COMPLEX*)gf, L, R, a, M, gdf);
    LTFAT_NAME_REAL(iwfacreal)((const LTFAT_COMPLEX*)gdf, L, R, a, M, gd);

#endif

error:
    LTFAT_SAFEFREEALL(gdf, gf);
    return status;
}


LTFAT_EXTERN int
LTFAT_NAME(gabtight_fir)(const LTFAT_TYPE* g, const ltfatInt gl,
                         const ltfatInt L, const ltfatInt a,
                         const ltfatInt M, const ltfatInt gtl, LTFAT_TYPE* gt)
{
    LTFAT_TYPE* tmpLong = NULL;

    int status = LTFATERR_SUCCESS;
    CHECKNULL(g); CHECKNULL(gt);
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive");
    CHECK(LTFATERR_NOTPOSARG, gtl > 0, "gtl must be positive");
    CHECK(LTFATERR_BADARG, L >= gl && L >= gtl,
          "L>=gl && L>= gtl must hold. Passed L=%d, gl=%d, gtl=%d", L, gl, gtl);

    CHECKMEM( tmpLong = ltfat_malloc(L * sizeof * tmpLong));

    CHECKSTATUS( LTFAT_NAME(fir2long)(g, gl, L, tmpLong), "fir2long failed");
    LTFAT_NAME(gabtight_long)(tmpLong, L, 1, a, M, tmpLong);
    CHECKSTATUS( LTFAT_NAME(long2fir)(tmpLong, L, gtl, gt), "long2fir failed");

error:
    if (tmpLong) ltfat_free(tmpLong);
    return status;
}
