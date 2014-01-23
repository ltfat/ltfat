#include "ltfat.h"
#include "ltfat_types.h"

#ifdef LTFAT_DOUBLE
#define LTFAT_POSV F77_FUNC(zposv,ZPOSV)
#define LTFAT_GESVD F77_FUNC(zgesvd,ZGESVD)
#define LTFAT_clapack_posv clapack_zposv
#define LTFAT_GEMM F77_FUNC (zgemm, ZGEMM)
#define LTFAT_cblas_gemm cblas_zgemm
#endif

#ifdef LTFAT_SINGLE
#define LTFAT_POSV F77_FUNC(cposv,CPOSV)
#define LTFAT_GESVD F77_FUNC(cgesvd,CGESVD)
#define LTFAT_clapack_posv clapack_cposv
#define LTFAT_GEMM F77_FUNC (cgemm, CGEMM)
#define LTFAT_cblas_gemm cblas_cgemm
#endif


/* Use the LAPACK library supplied with ATLAS */
#ifdef HAVE_ATLASLAPACK
#include "clapack.h"
#endif /* end of HAVE_ATLASLAPACK */

/* Call Fortran LAPACK */
#ifdef HAVE_LAPACK
#include "f77-fcn.h"


#ifdef __cplusplus
extern "C" {
#endif
F77_RET_T
LTFAT_POSV (F77_CONST_CHAR_ARG_DECL uplo,
            const ptrdiff_t *n, const ptrdiff_t *nrhs,
            const LTFAT_REAL *a, const ptrdiff_t *lda,
            LTFAT_REAL *b, const ptrdiff_t *ldb,
            ptrdiff_t *info
            F77_CHAR_ARG_LEN_DECL);


F77_RET_T
LTFAT_GESVD (F77_CONST_CHAR_ARG_DECL jobu,
             F77_CONST_CHAR_ARG_DECL jobvt,
             const ptrdiff_t *M, const ptrdiff_t *N,
             LTFAT_REAL* A, const ptrdiff_t *lda,
             LTFAT_REAL *S, LTFAT_REAL* U, const ptrdiff_t *ldu,
             LTFAT_REAL *VT, const ptrdiff_t *ldvt, LTFAT_REAL *work,
             const ptrdiff_t *lwork,
             LTFAT_REAL *rwork, ptrdiff_t *info
             F77_CHAR_ARG_LEN_DECL
             F77_CHAR_ARG_LEN_DECL);
#ifdef __cplusplus
}
#endif

#endif /* end of HAVE_LAPACK */

/* Call Fortran BLAS */
#ifdef HAVE_BLAS
#include "f77-fcn.h"

#ifdef __cplusplus
extern "C" {
#endif

F77_RET_T
LTFAT_GEMM (F77_CONST_CHAR_ARG_DECL TransA,
            F77_CONST_CHAR_ARG_DECL TransB,
            const ptrdiff_t *M, const ptrdiff_t *N,
            const ptrdiff_t *K,
            const LTFAT_COMPLEX *alpha,
            const LTFAT_COMPLEX *a, const ptrdiff_t *lda,
            const LTFAT_COMPLEX *b, const ptrdiff_t *ldb,
            const LTFAT_COMPLEX *beta,
            LTFAT_COMPLEX *c,
            const ptrdiff_t *ldc
            F77_CHAR_ARG_LEN_DECL
            F77_CHAR_ARG_LEN_DECL);
#ifdef __cplusplus
}
#endif

#endif /* end of HAVE_BLAS */



/* ----- Compute Cholesky factorization  ------------
 *
 * For simplification, the interface assumes that we
 * are using column-major format and storing the upper
 * triangle
 */
ltfatInt LTFAT_NAME(ltfat_posv)(const ptrdiff_t N, const ptrdiff_t NRHS,
                           LTFAT_COMPLEX *A, const ptrdiff_t lda,
                           LTFAT_COMPLEX *B, const ptrdiff_t ldb)
#ifdef HAVE_CBLASLAPACK
{
    return LTFAT_clapack_posv(CblasColMajor, CblasUpper, N, NRHS, (LTFAT_REAL*)A, lda,
                              (LTFAT_REAL*)B, ldb);
}
#endif
#ifdef HAVE_LAPACK
{
    ptrdiff_t info;
    char u;

    u = 'U';

    LTFAT_POSV (F77_CONST_CHAR_ARG2 (&u, 1),
                &N, &NRHS, (LTFAT_REAL*)A, &lda,
                (LTFAT_REAL*)B, &ldb,
                &info
                F77_CHAR_ARG_LEN (1)
               );

    return info;
}
#endif

/* ----- Compute SVD factorization  ------------ */
ltfatInt LTFAT_NAME(ltfat_gesvd)(const ptrdiff_t M, const ptrdiff_t N,
                            LTFAT_COMPLEX *A, const ptrdiff_t lda,
                            LTFAT_REAL *S, LTFAT_COMPLEX *U, const ptrdiff_t ldu,
                            LTFAT_COMPLEX *VT, const ptrdiff_t ldvt)
#ifdef HAVE_LAPACK
{

    ptrdiff_t lwork, info, maxMN;
    LTFAT_REAL workquery[2];
    LTFAT_REAL *rwork, *work;

    char jobu;
    char jobvt;

    /* Set constants to declare thin SVD */
    jobu = 'S';
    jobvt = 'S';

    maxMN = M > N ? M : N;

    /* Allocate workspace */
    rwork = (LTFAT_REAL*)ltfat_malloc(5*maxMN*sizeof(LTFAT_REAL));

    /* Ask ZGESVD what the dimension of WORK should be. */
    lwork = -1;
    LTFAT_GESVD (F77_CONST_CHAR_ARG2 (&jobu, 1),
                 F77_CONST_CHAR_ARG2 (&jobvt, 1),
                 &M, &N, (LTFAT_REAL*)A, &lda, S, (LTFAT_REAL*)U,
                 &ldu, (LTFAT_REAL*)VT,
                 &ldvt, (LTFAT_REAL*)(&workquery), &lwork,
                 rwork, &info
                 F77_CHAR_ARG_LEN (1)
                 F77_CHAR_ARG_LEN (1));

    /* Get the result from real part of work */
    lwork = (ptrdiff_t)(workquery[0]);

    /* Allocate more workspace */
    work = (LTFAT_REAL*)ltfat_malloc(lwork*sizeof(LTFAT_COMPLEX));

    /* Call the function for real this time */

    LTFAT_GESVD (F77_CONST_CHAR_ARG2 (&jobu, 1),
                 F77_CONST_CHAR_ARG2 (&jobvt, 1),
                 &M, &N, (LTFAT_REAL*)A, &lda, S,
                 (LTFAT_REAL*)U, &ldu, (LTFAT_REAL*)VT,
                 &ldvt, work, &lwork,
                 rwork, &info
                 F77_CHAR_ARG_LEN (1)
                 F77_CHAR_ARG_LEN (1));


    /* Free workspace */
    ltfat_free(rwork);
    ltfat_free(work);

    return info;

}
#endif /* End of HAVE_LAPACK */


void LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
                            const enum CBLAS_TRANSPOSE TransB,
                            const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
                            const LTFAT_COMPLEX *alpha,
                            const LTFAT_COMPLEX *A, const ptrdiff_t lda,
                            const LTFAT_COMPLEX *B, const ptrdiff_t ldb,
                            const LTFAT_COMPLEX *beta,
                            LTFAT_COMPLEX *C, const ptrdiff_t ldc)
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
#endif /* end of HAVE_BLAS */
