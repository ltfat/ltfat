#ifndef _GABDUALFIR_H
#define _GABDUALFIR_H

typedef enum { HANN, SQRTHANN, COS, SINE, HAMMING } LTFAT_FIRWIN;

#endif

LTFAT_EXTERN int
LTFAT_NAME(firwin)(LTFAT_FIRWIN win, int gl, LTFAT_TYPE* g);

LTFAT_EXTERN void
LTFAT_NAME(gabframediag)(const LTFAT_TYPE* g, ltfatInt gl,
                         ltfatInt a, ltfatInt M, ltfatInt dl, LTFAT_TYPE* d);

LTFAT_EXTERN int
LTFAT_NAME(gabtight_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                        ltfatInt M, LTFAT_TYPE* gt);

LTFAT_EXTERN int
LTFAT_NAME(gabdual_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a, ltfatInt M,
                       LTFAT_TYPE* gd);
