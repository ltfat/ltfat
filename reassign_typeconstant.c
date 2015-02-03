#include "ltfat.h"

LTFAT_EXTERN_TOO fbreassOptOut*
fbreassOptOut_init(const ltfatInt l, const ltfatInt inital)
{
   fbreassOptOut* ret = ltfat_calloc( 1, sizeof * ret);
   ret->l = l;
   // This is an array of pointers.
   ret->repos = ltfat_malloc(l * sizeof * ret->repos);
   ret->reposl = ltfat_calloc(l , sizeof * ret->reposl);
   ret->reposlmax = ltfat_malloc(l * sizeof * ret->reposlmax);
   ltfatInt inital2 = imax(1, inital);
   for (ltfatInt ii = 0; ii < l; ii++)
   {
      ret->repos[ii] = ltfat_malloc( inital2 * sizeof * ret->repos[ii]);
      ret->reposlmax[ii] = inital2;
   }

   return ret;
}

LTFAT_EXTERN_TOO void
fbreassOptOut_destroy(fbreassOptOut* oo)
{

   for (ltfatInt ii = 0; ii < oo->l; ii++)
   {
      if (oo->repos[ii] && oo->reposlmax[ii] > 0)
      {
         ltfat_free(oo->repos[ii]);
      }
   }

   LTFAT_SAFEFREEALL(oo->repos, oo->reposl, oo->reposlmax, oo);
   oo = NULL;
}

LTFAT_EXTERN_TOO void
fbreassOptOut_expand(fbreassOptOut* oo, ltfatInt ii)
{
   ltfatInt explmax = (ltfatInt) (fbreassOptOut_EXPANDRAT * oo->reposlmax[ii]);
   oo->repos[ii] = ltfat_realloc_and_copy(
                      oo->repos[ii],
                      oo->reposlmax[ii] * sizeof * oo->repos[ii],
                      explmax * sizeof * oo->repos[ii]);
   oo->reposlmax[ii] = explmax;
}

