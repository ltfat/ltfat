#ifndef _WAVELETS_H
#define _WAVELETS_H

LTFAT_EXTERN
ltfatExtType ltfatExtStringToEnum(const char* extType);


#endif

// CAN BE INCLUDED MORE THAN ONCE


LTFAT_EXTERN void
LTFAT_NAME(extend_left)(const LTFAT_TYPE *in,ltfatInt inLen, LTFAT_TYPE *buffer, ltfatInt buffLen, ltfatInt filtLen, ltfatExtType ext, ltfatInt a);

LTFAT_EXTERN void
LTFAT_NAME(extend_right)(const LTFAT_TYPE *in,ltfatInt inLen, LTFAT_TYPE *buffer, ltfatInt filtLen, ltfatExtType ext, ltfatInt a);




LTFAT_EXTERN void
LTFAT_NAME(convsub_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                       const ltfatInt L, const ltfatInt gl, const ltfatInt a, const ltfatInt skip,
                       LTFAT_TYPE *c, ltfatExtType ext);


LTFAT_EXTERN void
LTFAT_NAME(upconv_td)(const LTFAT_TYPE *c, const LTFAT_TYPE *g,
                      const ltfatInt L,  const ltfatInt gl, const ltfatInt a, const ltfatInt skip,
                      LTFAT_TYPE *f, ltfatExtType ext);


LTFAT_EXTERN void
LTFAT_NAME(filterbank_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g[],
                          const ltfatInt L, const ltfatInt gl[], const ltfatInt W,
                          const ltfatInt a[], const ltfatInt skip[], const ltfatInt M,
                          LTFAT_TYPE *c[], ltfatExtType ext);


LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_td)(const LTFAT_TYPE *c[], const LTFAT_TYPE *g[],
                           const ltfatInt L, const ltfatInt gl[], const ltfatInt W, const ltfatInt a[],
                           const ltfatInt skip[], const ltfatInt M, LTFAT_TYPE *f,
                           ltfatExtType ext);

LTFAT_EXTERN void
LTFAT_NAME(atrousfilterbank_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g[],
                                const ltfatInt L, const ltfatInt gl[], const ltfatInt W,
                                const ltfatInt a[], const ltfatInt skip[], const ltfatInt M,
                                LTFAT_TYPE *c, ltfatExtType ext);

LTFAT_EXTERN void
LTFAT_NAME(iatrousfilterbank_td)(const LTFAT_TYPE *c, const LTFAT_TYPE *g[],
                                 const ltfatInt L, const ltfatInt gl[], const ltfatInt W, const ltfatInt a[],
                                 const ltfatInt skip[], const ltfatInt M, LTFAT_TYPE *f,
                                 ltfatExtType ext);


LTFAT_EXTERN void
LTFAT_NAME(atrousconvsub_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                             const ltfatInt L, const ltfatInt gl,
                             const ltfatInt ga,ltfatInt skip,
                             LTFAT_TYPE *c, ltfatExtType ext);

LTFAT_EXTERN void
LTFAT_NAME(atrousupconv_td)(const LTFAT_TYPE *c, const LTFAT_TYPE *g,
                            const ltfatInt L, const ltfatInt gl,
                            const ltfatInt ga, const ltfatInt skip,
                            LTFAT_TYPE *f, ltfatExtType ext);








