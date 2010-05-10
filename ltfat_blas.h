#include "config.h"
#include "cblas.h"

void LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
			    const enum CBLAS_TRANSPOSE TransB,
			    const int M, const int N, const int K,
			    const LTFAT_COMPLEX *alpha,
			    const LTFAT_COMPLEX *A, const int lda,
			    const LTFAT_COMPLEX *B, const int ldb,
			    const LTFAT_COMPLEX *beta,
			    LTFAT_COMPLEX *C, const int ldc);

