#ifndef _LTFAT_REASSIGN_TYPECONSTANT_H
#define _LTFAT_REASSIGN_TYPECONSTANT_H

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

LTFAT_API fbreassOptOut*
fbreassOptOut_init(const ltfatInt l,const ltfatInt inital);

LTFAT_API void
fbreassOptOut_expand(fbreassOptOut* oo,const ltfatInt ii);

LTFAT_API void
fbreassOptOut_destroy(fbreassOptOut* oo);


#endif /* end of include guard: _REASSIGN_H */
