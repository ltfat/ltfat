#ifndef _RTDGT_H
#define _RTDGT_H

typedef enum
{
    LTFAT_RTDGTPHASE_ZERO,
    LTFAT_RTDGTPHASE_HALFSHIFT
} rtdgt_phasetype;

#endif

typedef struct
{
    const LTFAT_REAL* g;
    const ltfatInt gl;
    const ltfatInt M;
    const rtdgt_phasetype ptype;
    LTFAT_REAL* fftBuf;
    const fftBufLen;
    LTFAT_FFTW(plan) pfft;
} LTFAT_NAME(rtdgtreal_plan);

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                           const ltfatInt M, const rtdgt_phasetype ptype);

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_execute)(const LTFAT_NAME(rtdgtreal_plan)* p,
                              const LTFAT_REAL* f, const ltfatInt W,
                              LTFAT_COMPLEX* cout);

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_done)(const LTFAT_NAME(rtdgtreal_plan)* p);


LTFAT_EXTERN LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(irtdgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                            const ltfatInt M, const rtdgt_phasetype ptype);

LTFAT_EXTERN void
LTFAT_NAME(irtdgtreal_execute)(const LTFAT_NAME(rtdgtreal_plan)* p,
                               const LTFAT_COMPLEX* c, const ltfatInt W,
                               LTFAT_REAL* f);
