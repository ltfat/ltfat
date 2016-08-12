#include "ltfat/version.h"

static ltfat_library_version ltfat_version = { "0.1.0", __DATE__  " "  __TIME__,
                                               0, 1 , 0
                                             };

LTFAT_EXTERN ltfat_library_version* 
ltfat_get_library_version()
{
    return &ltfat_version;
}
