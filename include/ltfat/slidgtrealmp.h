#ifndef _LTFAT_SLIDGTREALMP_H
#define _LTFAT_SLIDGTREALMP_H


#endif

typedef struct LTFAT_NAME(slidgtrealmp_state) LTFAT_NAME(slidgtrealmp_state);

typedef int LTFAT_NAME(slidgtrealmp_niter_callback)(void* userdata,
        const LTFAT_REAL in[], int winLen, int taperLen, int zpadLen, int W, LTFAT_REAL out[]);

// LTFAT_API int
// LTFAT_NAME(slidgtrealmp_setnitercallback)(
//         LTFAT_NAME(slidgtrealmp_state)* p,
//         LTFAT_NAME(slidgtrealmp_niter_callback)* callback,
//         void* userdata);
//
//

LTFAT_API ltfat_int
LTFAT_NAME(slidgtrealmp_getprocdelay)( LTFAT_NAME(slidgtrealmp_state)* p);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_init)(
    LTFAT_NAME(dgtrealmp_parbuf)* pb, ltfat_int L,
    ltfat_int numChans, ltfat_int bufLenMax,
    LTFAT_NAME(slidgtrealmp_state)** pout);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_init_fromstates)(
    LTFAT_NAME(dgtrealmp_state)* mpstate,
    LTFAT_NAME(slicing_processor_state)* slistate,
    LTFAT_NAME(slidgtrealmp_state)** pout);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute)(
    LTFAT_NAME(slidgtrealmp_state)* p,
    const LTFAT_REAL* in[], ltfat_int inLen, ltfat_int chanNo,
    LTFAT_REAL* fout[]);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_execute_compact)(
    LTFAT_NAME(slidgtrealmp_state)* p,
    const LTFAT_REAL in[], ltfat_int inLen, ltfat_int chanNo,
    LTFAT_REAL fout[]);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_done)(LTFAT_NAME(slidgtrealmp_state)** p);

LTFAT_API int
LTFAT_NAME(slidgtrealmp_reset)(
        LTFAT_NAME(slidgtrealmp_state)* p);
