/*
Inplace friendly functions.
in==out is possibe.
*/
LTFAT_EXTERN void
LTFAT_NAME(dgtphaselockhelper)(LTFAT_TYPE *cin, const ltfatInt L,
                               const ltfatInt W, const ltfatInt a,
                               const ltfatInt M, LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtphaseunlockhelper)(LTFAT_TYPE *cin, const ltfatInt L,
                                 const ltfatInt W, const ltfatInt a,
                                 const ltfatInt M, LTFAT_TYPE *cout);
LTFAT_EXTERN
void LTFAT_NAME(circshift)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ltfatInt L, const ltfatInt shift);

LTFAT_EXTERN
void LTFAT_NAME(reverse_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ltfatInt L);

LTFAT_EXTERN
void LTFAT_NAME(conjugate_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ltfatInt L);

LTFAT_EXTERN
void LTFAT_NAME(periodize_array)(LTFAT_TYPE *in, const ltfatInt Lin,
                                 LTFAT_TYPE *out, const ltfatInt Lout);

LTFAT_EXTERN
void LTFAT_NAME(findmaxinarray)(const LTFAT_TYPE *in, const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx);

LTFAT_EXTERN
int LTFAT_NAME(findmaxinarraywrtmask)(const LTFAT_TYPE *in, const int *mask, const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx);

LTFAT_EXTERN
void LTFAT_NAME(array2complex)(LTFAT_TYPE *in, LTFAT_COMPLEX *out, const ltfatInt L);
