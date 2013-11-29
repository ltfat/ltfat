#include "wavelets.h"
#include "goertzel.h"
#include "ciutils.h"

LTFAT_EXTERN void
LTFAT_NAME(col2diag)(const LTFAT_TYPE *cin, const int L,
		       LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(gabdual_long)(const LTFAT_TYPE *g,
				const int L, const int R, const int a,
				const int M, LTFAT_TYPE *gd);

LTFAT_EXTERN void
LTFAT_NAME(gabtight_long)(const LTFAT_TYPE *g,
				 const int L, const int R, const int a,
				 const int M, LTFAT_TYPE *gd);


/* --------- Wilson and WMDCT bases ---------*/
LTFAT_EXTERN void
LTFAT_NAME(dwilt_long)(const LTFAT_TYPE *f,
			     const LTFAT_TYPE *g,
			     const int L, const int W, const int M,
			     LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dwilt_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
			    const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_TYPE *cout);


LTFAT_EXTERN void
LTFAT_NAME(dwiltiii_long)(const LTFAT_TYPE *f,
			     const LTFAT_TYPE *g,
			     const Lint L, const Lint W, const Lint M,
			     LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dwiltiii_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
			    const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_TYPE *cout);


/* --------- Wilson and WMDCT inverses ---------*/


LTFAT_EXTERN void
LTFAT_NAME(idwilt_long)(const LTFAT_TYPE *cin,
			     const LTFAT_TYPE *g,
			     const int L, const int W, const int M,
			     LTFAT_TYPE *f);

LTFAT_EXTERN void
LTFAT_NAME(idwilt_fb)(const LTFAT_TYPE *cin, const LTFAT_TYPE *g,
			    const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_TYPE *f);

LTFAT_EXTERN void
LTFAT_NAME(idwiltiii_long)(const LTFAT_TYPE *cin,
			     const LTFAT_TYPE *g,
			     const Lint L, const Lint W, const Lint M,
			     LTFAT_TYPE *f);

LTFAT_EXTERN void
LTFAT_NAME(idwiltiii_fb)(const LTFAT_TYPE *cin, const LTFAT_TYPE *g,
			    const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_TYPE *f);


