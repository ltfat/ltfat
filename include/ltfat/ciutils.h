#ifndef _ciutils_h
#define _ciutils_h

typedef enum
{
    LTFAT_NORMALIZE_NULL = 0,
    LTFAT_NORMALIZE_AREA = 1,
    LTFAT_NORMALIZE_ENERGY = 2,
    LTFAT_NORMALIZE_RMS = 3,
    LTFAT_NORMALIZE_WAV = 4,
    LTFAT_NORMALIZE_INF
} ltfat_normalize_t;

#endif


LTFAT_EXTERN void
LTFAT_NAME(dgtphaselockhelper)(LTFAT_TYPE *cin, const ltfatInt L,
                               const ltfatInt W, const ltfatInt a,
                               const ltfatInt M, LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtphaseunlockhelper)(LTFAT_TYPE *cin, const ltfatInt L,
                                 const ltfatInt W, const ltfatInt a,
                                 const ltfatInt M, LTFAT_TYPE *cout);
LTFAT_EXTERN
void LTFAT_NAME(circshift)(const LTFAT_TYPE *in, const ltfatInt L,
                           const ltfatInt shift, LTFAT_TYPE *out);

LTFAT_EXTERN void
LTFAT_NAME(fftshift)(const LTFAT_TYPE* in, ltfatInt L, LTFAT_TYPE* out);

LTFAT_EXTERN void
LTFAT_NAME(ifftshift)(const LTFAT_TYPE* in, ltfatInt L, LTFAT_TYPE* out);

LTFAT_EXTERN
void LTFAT_NAME(reverse_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ltfatInt L);

LTFAT_EXTERN
void LTFAT_NAME(conjugate_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ltfatInt L);

LTFAT_EXTERN
void LTFAT_NAME(periodize_array)(LTFAT_TYPE *in, const ltfatInt Lin,
                                 LTFAT_TYPE *out, const ltfatInt Lout);

LTFAT_EXTERN
void LTFAT_NAME(fold_array)(const LTFAT_TYPE *in, const ltfatInt Lin,
                            const ltfatInt Lfold, LTFAT_TYPE *out);

LTFAT_EXTERN
void LTFAT_NAME(findmaxinarray)(const LTFAT_TYPE *in, const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx);

LTFAT_EXTERN int
LTFAT_NAME(findmaxinarraywrtmask)(const LTFAT_TYPE *in, const int *mask, const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx);

LTFAT_EXTERN
void LTFAT_NAME(array2complex)(LTFAT_TYPE *in, LTFAT_COMPLEX *out, const ltfatInt L);

LTFAT_EXTERN void
LTFAT_NAME(fir2long)(const LTFAT_TYPE *f, const ltfatInt Lfir, const ltfatInt Llong,
                     LTFAT_TYPE *h);

LTFAT_EXTERN void
LTFAT_NAME(long2fir)(const LTFAT_TYPE *f, const ltfatInt Llong, const ltfatInt Lfir,
                     LTFAT_TYPE *h);
