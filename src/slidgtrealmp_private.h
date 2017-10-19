#ifndef _LTFAT_SLIDGTREALMP_PRIVATE_H
#define _LTFAT_SLIDGTREALMP_PRIVATE_H
#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealmp_private.h"

struct LTFAT_NAME(slidgtrealmp_state)
{
    LTFAT_NAME(dgtrealmp_state)* istate;
    int owningistate;
};



#endif
