#ifndef _LTFAT_SLIDGTREALMP_PRIVATE_H
#define _LTFAT_SLIDGTREALMP_PRIVATE_H


#endif

struct LTFAT_NAME(slidgtrealmp_state)
{
    LTFAT_NAME(dgtrealmp_state)* mpstate;
    LTFAT_NAME(slicing_processor_state)* slistate;
    LTFAT_COMPLEX** couttmp;
};

int
LTFAT_NAME(slidgtrealmp_execute_callback)(void* userdata,
        const LTFAT_REAL in[], int winLen, int taperLen, 
        int zpadLen, int W, LTFAT_REAL out[]);
