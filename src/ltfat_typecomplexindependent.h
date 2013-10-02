#include "wavelets.h"

LTFAT_EXTERN void
LTFAT_H_NAME(col2diag)(const LTFAT_H_TYPE *cin, const int L,
		       LTFAT_H_TYPE *cout);

LTFAT_EXTERN
void LTFAT_H_NAME(gabdual_long)(const LTFAT_H_TYPE *g,
				const int L, const int R, const int a,
				const int M, LTFAT_H_TYPE *gd);

LTFAT_EXTERN
void LTFAT_H_NAME(gabtight_long)(const LTFAT_H_TYPE *g,
				 const int L, const int R, const int a,
				 const int M, LTFAT_H_TYPE *gd);


LTFAT_EXTERN
void LTFAT_H_NAME(gga)(const LTFAT_H_TYPE *fPtr, const double *indVecPtr,
                  const int L, const int W, const int M,
		          LTFAT_H_COMPLEXH *cPtr);

LTFAT_EXTERN
void LTFAT_H_NAME(circshift)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const ptrdiff_t L, const ptrdiff_t shift);

LTFAT_EXTERN
void LTFAT_H_NAME(reverse_array)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const size_t L);





