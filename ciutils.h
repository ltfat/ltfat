/*
Inplace friendly functions.
in==out is possibe.
*/

LTFAT_EXTERN
void LTFAT_NAME(circshift)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ptrdiff_t L, const ptrdiff_t shift);

LTFAT_EXTERN
void LTFAT_NAME(reverse_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, const size_t L);

LTFAT_EXTERN
void LTFAT_NAME(conjugate_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, const size_t L);

LTFAT_EXTERN
void LTFAT_NAME(array2complex)(LTFAT_TYPE *in, LTFAT_COMPLEX *out, const size_t L);
