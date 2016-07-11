#ifndef _LTFAT_VERSION_H
#define _LTFAT_VERSION_H


typedef struct
{
    const char* version;
    const char* build_date;
} ltfat_library_version;


ltfat_library_version*
get_ltfat_library_version();


#endif
