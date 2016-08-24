#include "ltfat/version.h"

#define LTFAT_MAKEVESRIONSTRING(major,minor,patch) #major "." #minor "." #patch
#define LTFAT_VERSIONSTRING(major,minor,patch) LTFAT_MAKEVESRIONSTRING(major,minor,patch)


static ltfat_library_version ltfat_version =
{
    LTFAT_VERSIONSTRING(LTFAT_VERSION_MAJOR,LTFAT_VERSION_MINOR,LTFAT_VERSION_PATCH),
    __DATE__  " "  __TIME__,
    LTFAT_VERSION_MAJOR, LTFAT_VERSION_MINOR , LTFAT_VERSION_PATCH
};

LTFAT_API ltfat_library_version*
ltfat_get_version()
{
    return &ltfat_version;
}

LTFAT_API int
ltfat_is_compatible_version()
{
    ltfat_library_version* dll_version = ltfat_get_version();
    return (int) ( dll_version->major == LTFAT_VERSION_MAJOR);
}
