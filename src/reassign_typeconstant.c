#include "ltfat.h"
#include "ltfat/macros.h"

LTFAT_API fbreassOptOut*
fbreassOptOut_init(const ltfatInt l, const ltfatInt inital)
{
    fbreassOptOut* ret = (fbreassOptOut*) ltfat_calloc( 1, sizeof * ret);
    ret->l = l;
    // This is an array of pointers.
    ret->repos = (ltfatInt**) ltfat_malloc(l * sizeof * ret->repos);
    ret->reposl = (ltfatInt*) ltfat_calloc(l , sizeof * ret->reposl);
    ret->reposlmax = (ltfatInt*) ltfat_malloc(l * sizeof * ret->reposlmax);
    ltfatInt inital2 = ltfat_imax(1, inital);
    for (ltfatInt ii = 0; ii < l; ii++)
    {
        ret->repos[ii] = (ltfatInt*) ltfat_malloc( inital2 * sizeof * ret->repos[ii]);
        ret->reposlmax[ii] = inital2;
    }

    return ret;
}

LTFAT_API void
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

LTFAT_API void
fbreassOptOut_expand(fbreassOptOut* oo, ltfatInt ii)
{
    ltfatInt explmax = (ltfatInt) (fbreassOptOut_EXPANDRAT * oo->reposlmax[ii]);
    oo->repos[ii] = (ltfatInt*) ltfat_realloc( (void*) oo->repos[ii],
                    oo->reposlmax[ii] * sizeof * oo->repos[ii],
                    explmax * sizeof * oo->repos[ii]);
    oo->reposlmax[ii] = explmax;
}
