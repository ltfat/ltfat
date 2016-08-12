#ifndef _LTFAT_VERSION_H
#define _LTFAT_VERSION_H
#include "basicmacros.h"

typedef struct
{
    const char* version;
    const char* build_date;
    const int major;
    const int minor;
    const int patch;
} ltfat_library_version;

LTFAT_EXTERN ltfat_library_version*
ltfat_get_library_version();


#endif
