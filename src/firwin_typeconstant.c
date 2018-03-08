#include "ltfat.h"

LTFAT_API int
ltfat_str2firwin(const char* win)
{
    if( strcmp("hann",win) || strcmp("hanning",win) || strcmp("nuttall10",win))
        return LTFAT_HANN;
    else if( strcmp("sqrthann",win) || strcmp("cosine",win) || strcmp("sine",win))
        return LTFAT_SQRTHANN;
    else if( strcmp("nuttall01",win)) return LTFAT_NUTTALL01;
    else if( strcmp("square",win) || strcmp("rect",win)) return LTFAT_SQUARE;
    else if( strcmp("tria",win) || strcmp("triangular",win) || strcmp("bartlett",win))
        return LTFAT_TRIA;
    else if( strcmp("sqrttria",win)) return LTFAT_SQRTTRIA;
    else if( strcmp("blackman",win)) return LTFAT_BLACKMAN;
    else if( strcmp("blackman2",win)) return LTFAT_BLACKMAN2;
    else if( strcmp("nuttall",win) || strcmp("nuttall12",win)) return LTFAT_NUTTALL;
    else if( strcmp("ogg",win) || strcmp("itersine",win)) return LTFAT_OGG;
    else if( strcmp("nuttall20",win)) return LTFAT_NUTTALL20;
    else if( strcmp("nuttall11",win)) return LTFAT_NUTTALL11;
    else if( strcmp("nuttall02",win)) return LTFAT_NUTTALL02;
    else if( strcmp("nuttall30",win)) return LTFAT_NUTTALL30;
    else if( strcmp("nuttall21",win)) return LTFAT_NUTTALL21;
    else if( strcmp("nuttall03",win)) return LTFAT_NUTTALL03;
    else if( strcmp("truncgauss01",win)) return LTFAT_TRUNCGAUSS01;

    return LTFATERR_BADARG;
}
