/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE

#include "config.h"
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(gabtight_long)(const LTFAT_TYPE *g,
			  const int L, const int R, const int a,
			  const int M, LTFAT_TYPE *gd)
{

#ifdef LTFAT_COMPLEXTYPE

   LTFAT_COMPLEX *gf = ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *gdf = ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));

   LTFAT_NAME_REAL(wfac)(g, L, R, a, M, gf);
   LTFAT_NAME_REAL(gabtight_fac)((const LTFAT_COMPLEX *)gf,L,R,a,M,gdf);
   LTFAT_NAME_REAL(iwfac)((const LTFAT_COMPLEX *)gdf,L,R,a,M,gd);

#else

   const int wfs = L; /* wfacreal_size(L,a,M); */

   LTFAT_COMPLEX *gf = ltfat_malloc(wfs*R*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *gdf = ltfat_malloc(wfs*R*sizeof(LTFAT_COMPLEX));

   LTFAT_NAME_REAL(wfacreal)(g, L, R, a, M, gf);
   LTFAT_NAME_REAL(gabtightreal_fac)((const LTFAT_COMPLEX *)gf,L,R,a,M,gdf);
   LTFAT_NAME_REAL(iwfacreal)((const LTFAT_COMPLEX *)gdf,L,R,a,M,gd);

#endif

   ltfat_free(gdf);
   ltfat_free(gf);


}
#endif
