/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE

#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(gabdual_long)(const LTFAT_TYPE *g,
                         const ltfatInt L, const ltfatInt R, const ltfatInt a,
                         const ltfatInt M, LTFAT_TYPE *gd)
{

    const ltfatInt wfs = L;


    LTFAT_COMPLEX *gf = ltfat_malloc(wfs*R*sizeof(LTFAT_COMPLEX));
    LTFAT_COMPLEX *gdf = ltfat_malloc(wfs*R*sizeof(LTFAT_COMPLEX));

#ifdef LTFAT_COMPLEXTYPE

    LTFAT_NAME(wfac)(g, L, R, a, M, gf);
    LTFAT_NAME_REAL(gabdual_fac)((const LTFAT_COMPLEX *)gf,L,R,a,M,gdf);
    LTFAT_NAME(iwfac)((const LTFAT_COMPLEX *)gdf,L,R,a,M,gd);

#else

    LTFAT_NAME_REAL(wfacreal)(g, L, R, a, M, gf);
    LTFAT_NAME_REAL(gabdualreal_fac)((const LTFAT_COMPLEX *)gf,L,R,a,M,gdf);
    LTFAT_NAME_REAL(iwfacreal)((const LTFAT_COMPLEX *)gdf,L,R,a,M,gd);

#endif

    LTFAT_SAFEFREEALL(gdf,gf);
}

#endif
