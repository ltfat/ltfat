#ifndef _LTFAT_VERSION_H
#define _LTFAT_VERSION_H
#include "basicmacros.h"

#define LTFAT_VERSION_MAJOR 0
#define LTFAT_VERSION_MINOR 1
#define LTFAT_VERSION_PATCH 0

typedef struct
{
    const char* version;
    const char* build_date;
    const int major;
    const int minor;
    const int patch;
} ltfat_library_version;


#ifdef __cplusplus
extern "C"
{
#endif

LTFAT_API ltfat_library_version*
ltfat_get_version();

LTFAT_API int
ltfat_is_compatible_version();

#ifdef __cplusplus
}  // extern "C"
#endif


#endif
