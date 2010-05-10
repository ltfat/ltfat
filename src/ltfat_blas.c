#include "config.h"
#include "cblas.h"

#ifdef LTFAT_DOUBLE
#define LTFAT_GEMM F77_FUNC (zgemm, ZGEMM)
#define LTFAT_cblas_gemm cblas_zgemm
#endif

#ifdef LTFAT_SINGLE
#define LTFAT_GEMM F77_FUNC (cgemm, CGEMM)
#define LTFAT_cblas_gemm cblas_cgemm
#endif

#ifdef HAVE_BLAS
#include "f77-fcn.h"

#ifdef __cplusplus
extern "C" {
#endif

  F77_RET_T
  LTFAT_GEMM (F77_CONST_CHAR_ARG_DECL TransA,
	      F77_CONST_CHAR_ARG_DECL TransB,
	      const int *M, const int *N,
	      const int *K, 
	      const LTFAT_COMPLEX *alpha,
	      const LTFAT_COMPLEX *a, const int *lda,
	      const LTFAT_COMPLEX *b, const int *ldb, 
	      const LTFAT_COMPLEX *beta,
	      LTFAT_COMPLEX *c,
	      const int *ldc
	      F77_CHAR_ARG_LEN_DECL
	      F77_CHAR_ARG_LEN_DECL);
#ifdef __cplusplus
}
#endif

#endif /* end of HAVE_BLAS */


void LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
		 const int M, const int N, const int K,
		 const LTFAT_COMPLEX *alpha,
		 const LTFAT_COMPLEX *A, const int lda,
		 const LTFAT_COMPLEX *B, const int ldb,
                 const LTFAT_COMPLEX *beta,
		 LTFAT_COMPLEX *C, const int ldc)
#ifdef HAVE_CBLAS
{

  LTFAT_cblas_gemm(CblasColMajor, TransA, TransB, M, N, K,
	      (LTFAT_REAL*)alpha, (LTFAT_REAL*)A, lda, (LTFAT_REAL*)B, ldb,
	      (LTFAT_REAL*)beta, (LTFAT_REAL*)C, ldc);

}
#endif
#ifdef HAVE_BLAS
{
  char ca, cb;

  if (TransA == CblasNoTrans)   ca='N';
  if (TransA == CblasConjTrans) ca='C';

  if (TransB == CblasNoTrans)   cb='N';
  if (TransB == CblasConjTrans) cb='C';

  LTFAT_GEMM (F77_CONST_CHAR_ARG2 (&ca, 1),
	      F77_CONST_CHAR_ARG2 (&cb, 1),
	      &M, &N, &K,
	      alpha,
	      A, &lda,
	      B, &ldb,
	      beta, C, &ldc
	      F77_CHAR_ARG_LEN (1)
	      F77_CHAR_ARG_LEN (1)
	      );
  
  
}
#endif
