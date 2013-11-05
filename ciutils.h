/*
Inplace friendly functions.
in==out is possibe.
*/

LTFAT_EXTERN
void LTFAT_H_NAME(circshift)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const ptrdiff_t L, const ptrdiff_t shift);

LTFAT_EXTERN
void LTFAT_H_NAME(reverse_array)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const size_t L);

LTFAT_EXTERN
void LTFAT_H_NAME(conjugate_array)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const size_t L);

LTFAT_EXTERN
void LTFAT_H_NAME(array2complex)(LTFAT_H_TYPE *in, LTFAT_H_COMPLEXH *out, const size_t L);
