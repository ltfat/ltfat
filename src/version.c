#include "ltfat/version.h"

static ltfat_library_version ltfat_version = { "0.1.0", __DATE__  " "  __TIME__,
                                               0, 1 , 0
                                             };

ltfat_library_version* get_ltfat_library_version()
{
    return &ltfat_version;
}
