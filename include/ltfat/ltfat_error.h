#ifndef _LTFAT_ERROR_H
#define _LTFAT_ERROR_H

enum
{
    LTFATERR_SUCCESS     =  0,
    LTFATERR_FAILURE     = -1,
    LTFATERR_MEMERR      =  1,
    LTFATERR_UNKNOWNWIN  =  2,
    LTFATERR_NOTAFRAME   =  3,
    LTFATERR_NOTPAINLESS =  4
};


#ifdef __cplusplus
extern "C"
{
#endif

const char*
ltfat_strerror(const int err);
// {
//     switch (err)
//     {
//     case LTFAT_SUCCESS:
//         return "success";
//     default:
//         return "unknown error code";
//     };
//
// }

#ifdef __cplusplus
}
#endif



#endif
