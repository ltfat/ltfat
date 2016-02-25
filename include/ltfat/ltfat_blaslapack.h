#ifndef _ltfat_blaslapack
#define _ltfat_blaslapack
#include "cblas.h"

// LAPACK overwrites the input argument.
ltfatInt
LTFAT_NAME(ltfat_posv)(const ptrdiff_t N, const ptrdiff_t NRHS,
                       LTFAT_COMPLEX *A, const ptrdiff_t lda,
                       LTFAT_COMPLEX *B, const ptrdiff_t ldb);

// LAPACK overwrites the input argument.
ltfatInt
LTFAT_NAME(ltfat_gesvd)(const ptrdiff_t M, const ptrdiff_t N,
                        LTFAT_COMPLEX *A, const ptrdiff_t lda,
                        LTFAT_REAL *S, LTFAT_COMPLEX *U, const ptrdiff_t ldu,
                        LTFAT_COMPLEX *VT, const ptrdiff_t ldvt);

void
LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
                       const enum CBLAS_TRANSPOSE TransB,
                       const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
                       const LTFAT_COMPLEX *alpha,
                       const LTFAT_COMPLEX *A, const ptrdiff_t lda,
                       const LTFAT_COMPLEX *B, const ptrdiff_t ldb,
                       const LTFAT_COMPLEX *beta,
                       LTFAT_COMPLEX *C, const ptrdiff_t ldc);

#endif
