#include "wavelets.h"
#include "goertzel.h"
#include "ciutils.h"

LTFAT_EXTERN void
LTFAT_NAME(col2diag)(const LTFAT_TYPE *cin, const ltfatInt L,
                     LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(gabdual_long)(const LTFAT_TYPE *g,
                         const ltfatInt L, const ltfatInt R, const ltfatInt a,
                         const ltfatInt M, LTFAT_TYPE *gd);

LTFAT_EXTERN void
LTFAT_NAME(gabtight_long)(const LTFAT_TYPE *g,
                          const ltfatInt L, const ltfatInt R, const ltfatInt a,
                          const ltfatInt M, LTFAT_TYPE *gd);


/* --------- Wilson and WMDCT bases ---------*/
LTFAT_EXTERN void
LTFAT_NAME(dwilt_long)(const LTFAT_TYPE *f,
                       const LTFAT_TYPE *g,
                       const ltfatInt L, const ltfatInt W, const ltfatInt M,
                       LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dwilt_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                     const ltfatInt L, const ltfatInt gl, const ltfatInt W, const ltfatInt M,
                     LTFAT_TYPE *cout);


LTFAT_EXTERN void
LTFAT_NAME(dwiltiii_long)(const LTFAT_TYPE *f,
                          const LTFAT_TYPE *g,
                          const ltfatInt L, const ltfatInt W, const ltfatInt M,
                          LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dwiltiii_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                        const ltfatInt L, const ltfatInt gl, const ltfatInt W, const ltfatInt M,
                        LTFAT_TYPE *cout);


/* --------- Wilson and WMDCT inverses ---------*/


LTFAT_EXTERN void
LTFAT_NAME(idwilt_long)(const LTFAT_TYPE *cin,
                        const LTFAT_TYPE *g,
                        const ltfatInt L, const ltfatInt W, const ltfatInt M,
                        LTFAT_TYPE *f);

LTFAT_EXTERN void
LTFAT_NAME(idwilt_fb)(const LTFAT_TYPE *cin, const LTFAT_TYPE *g,
                      const ltfatInt L, const ltfatInt gl, const ltfatInt W, const ltfatInt M,
                      LTFAT_TYPE *f);

LTFAT_EXTERN void
LTFAT_NAME(idwiltiii_long)(const LTFAT_TYPE *cin,
                           const LTFAT_TYPE *g,
                           const ltfatInt L, const ltfatInt W, const ltfatInt M,
                           LTFAT_TYPE *f);

LTFAT_EXTERN void
LTFAT_NAME(idwiltiii_fb)(const LTFAT_TYPE *cin, const LTFAT_TYPE *g,
                         const ltfatInt L, const ltfatInt gl, const ltfatInt W, const ltfatInt M,
                         LTFAT_TYPE *f);

/* --------------- DCT -------------------*/

LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(dct_init)( const ltfatInt L, const ltfatInt W, LTFAT_TYPE *cout,
                      const dct_kind kind);


LTFAT_EXTERN void
LTFAT_NAME(dct)(const LTFAT_TYPE *f, const ltfatInt L, const ltfatInt W,
                LTFAT_TYPE *cout, const dct_kind kind);

LTFAT_EXTERN void
LTFAT_NAME(dct_execute)(const LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
                        const ltfatInt L, const ltfatInt W,
                        LTFAT_TYPE *cout, const dct_kind kind);

/* --------------- DST -------------------*/

LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(dst_init)( const ltfatInt L, const ltfatInt W, LTFAT_TYPE *cout,
                      const dst_kind kind);

LTFAT_EXTERN void
LTFAT_NAME(dst)(const LTFAT_TYPE *f, const ltfatInt L, const ltfatInt W,
                LTFAT_TYPE *cout, const dst_kind kind);

LTFAT_EXTERN void
LTFAT_NAME(dst_execute)(LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
                        const ltfatInt L, const ltfatInt W, LTFAT_TYPE *cout,
                        const dst_kind kind);


