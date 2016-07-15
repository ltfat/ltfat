#ifndef _LTFAT_VERSION_H
#define _LTFAT_VERSION_H


typedef struct
{
    const char* version;
    const char* build_date;
    const int major;
    const int minor;
    const int patch;
} ltfat_library_version;


ltfat_library_version*
get_ltfat_library_version();


#endif
