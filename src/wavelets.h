#ifndef _WAVELETS_H
#define _WAVELETS_H
#include <string.h>



inline static ltfatExtType ltfatExtStringToEnum(const char* extType)
{
    if(strcmp(extType,"per")==0)
    {
        return PER;
    }
    else if(strcmp(extType,"perdec")==0)
    {
        return PERDEC;
    }
    else if(strcmp(extType,"ppd")==0)
    {
        return PPD;
    }
    else if(strcmp(extType,"sym")==0)
    {
        return SYM;
    }
    else if(strcmp(extType,"even")==0)
    {
        return EVEN;
    }
    else if(strcmp(extType,"symw")==0)
    {
        return SYMW;
    }
    else if(strcmp(extType,"odd")==0)
    {
        return ODD;
    }
    else if(strcmp(extType,"asymw")==0)
    {
        return ASYMW;
    }
    else if(strcmp(extType,"sp0")==0)
    {
        return SP0;
    }
    else if(strcmp(extType,"zpd")==0)
    {
        return ZPD;
    }
    else if(strcmp(extType,"zero")==0)
    {
        return ZERO;
    }
    else if(strcmp(extType,"valid")==0)
    {
        return VALID;
    }
    else
    {
        return BAD_TYPE;
    }
}





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








