/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE

#include "config.h"
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(gabdual_long)(const LTFAT_TYPE *g,
			 const int L, const int R, const int a,
			 const int M, LTFAT_TYPE *gd)
{

#ifdef LTFAT_COMPLEXTYPE
   const int wfs = L;
#else
   const int wfs = L; /* wfacreal_size(L,a,M); */
#endif

   LTFAT_COMPLEX *gf = ltfat_malloc(wfs*R*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *gdf = ltfat_malloc(wfs*R*sizeof(LTFAT_COMPLEX));

#ifdef LTFAT_COMPLEXTYPE
  
   LTFAT_NAME_REAL(wfac)(g, L, R, a, M, gf);
   LTFAT_NAME_REAL(gabdual_fac)((const LTFAT_COMPLEX *)gf,L,R,a,M,gdf);
   LTFAT_NAME_REAL(iwfac)((const LTFAT_COMPLEX *)gdf,L,R,a,M,gd);

#else
    
   LTFAT_NAME_REAL(wfacreal)(g, L, R, a, M, gf);
   LTFAT_NAME_REAL(gabdualreal_fac)((const LTFAT_COMPLEX *)gf,L,R,a,M,gdf);
   LTFAT_NAME_REAL(iwfacreal)((const LTFAT_COMPLEX *)gdf,L,R,a,M,gd);

#endif

   ltfat_free(gdf);
   ltfat_free(gf);

}

#endif
