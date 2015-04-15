#ifndef _REASSIGN_TYPECONSTANT_H
#define _REASSIGN_TYPECONSTANT_H
#include "ltfat.h"

#define fbreassOptOut_EXPANDRAT 2

typedef enum
{
    REASS_DEFAULT          = 0,
    REASS_NOTIMEWRAPAROUND = 1
} fbreassHints;

typedef struct {
   ltfatInt** repos;
   ltfatInt*  reposl;
   ltfatInt*  reposlmax;
   ltfatInt   l;
} fbreassOptOut;

LTFAT_EXTERN_TOO fbreassOptOut*
fbreassOptOut_init(const ltfatInt l,const ltfatInt inital);

LTFAT_EXTERN_TOO void
fbreassOptOut_expand(fbreassOptOut* oo,const ltfatInt ii);

LTFAT_EXTERN_TOO void
fbreassOptOut_destroy(fbreassOptOut* oo);


#endif /* end of include guard: _REASSIGN_H */
