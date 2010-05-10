#include "config.h"
#include "cblas.h"

/* LAPACK overwrites the input argument. */
int LTFAT_NAME(ltfat_posv)(const int N, const int NRHS,
	       LTFAT_COMPLEX *A, const int lda,
	       LTFAT_COMPLEX *B, const int ldb);

/* LAPACK overwrites the input argument. */
int LTFAT_NAME(ltfat_gesvd)(const int M, const int N,
		 LTFAT_COMPLEX *A, const int lda,
		 LTFAT_REAL *S, LTFAT_COMPLEX *U, const int ldu,
		 LTFAT_COMPLEX *VT, const int ldvt);
