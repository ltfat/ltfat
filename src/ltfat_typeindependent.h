/*  --------- factorizations --------------- */

extern void
LTFAT_H_NAME(wfac)(const LTFAT_H_COMPLEX *g, const int L, 
		 const int a, const int M, LTFAT_H_COMPLEX *gf);
		
extern void
LTFAT_H_NAME(wfac_r)(const LTFAT_H_REAL *g, const int L, 
		   const int a, const int M, LTFAT_H_COMPLEX *gf);

extern void
LTFAT_H_NAME(iwfac)(const LTFAT_H_COMPLEX *gf, const int L, 
		  const int a, const int M, LTFAT_H_COMPLEX *g);

extern void
LTFAT_H_NAME(iwfac_r)(const LTFAT_H_COMPLEX *gf, const int L, 
		    const int a, const int M, LTFAT_H_REAL *g);

/* --------- DGT by factorization ------------ */

extern void
LTFAT_H_NAME(dgt_fac)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *gf,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(dgt_long)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(dgtreal_long)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_H_COMPLEX *cout);
		
extern void
LTFAT_H_NAME(dgt_fac_r)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf,
			const int L,
			const int W,  const int a,
			const int M, LTFAT_H_COMPLEX *cout);
							
extern
void LTFAT_H_NAME(dgtreal_fac)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf,
			       const int L,
			       const int W,  const int a,
			       const int M, LTFAT_H_COMPLEX *cout);
	
extern void
LTFAT_H_NAME(dgt_walnut)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *gf,
		       const int L, const int W,
		        const int a, const int M, LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(dgt_walnut_r)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf, 
			 const int L, const int W,
			  const int a, const int M, LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(dgtreal_walnut)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf, 
			   const int L, const int W,
			    const int a, const int M, LTFAT_H_REAL *cout);

extern void
LTFAT_H_NAME(idgt_fac)(const LTFAT_H_COMPLEX *c, const LTFAT_H_COMPLEX *gf,
		       const int L,
		       const int W,const int a, const int M,
		       LTFAT_H_COMPLEX *f);

extern void 
LTFAT_H_NAME(idgtreal_fac)(const LTFAT_H_COMPLEX *cin, const LTFAT_H_COMPLEX *gf, 
			   const int L, const int W,
			   const int a, const int M,
			   LTFAT_H_REAL *f);


/* --------- Wilson and WMDCT bases ---------*/
void LTFAT_H_NAME(dwilt_long)(const LTFAT_H_COMPLEX *f,
			     const LTFAT_H_COMPLEX *g,
			     const int L, const int W, const int M, 
			     LTFAT_H_COMPLEX *cout);

void LTFAT_H_NAME(dwiltreal_long)(const LTFAT_H_REAL *f,
			     const LTFAT_H_REAL *g,
			     const int L, const int W, const int M, 
			     LTFAT_H_REAL *cout);

void LTFAT_H_NAME(wmdct_long)(const LTFAT_H_COMPLEX *f,
			     const LTFAT_H_COMPLEX *g,
			     const int L, const int W, const int M, 
			     LTFAT_H_COMPLEX *cout);

void LTFAT_H_NAME(wmdctreal_long)(const LTFAT_H_REAL *f,
			     const LTFAT_H_REAL *g,
			     const int L, const int W, const int M, 
			     LTFAT_H_REAL *cout);

/* --------- dual windows etc. --------------- */

extern void
LTFAT_H_NAME(gabdual_fac)(const LTFAT_H_COMPLEX *g, const int L, 
			const int a, const int M, LTFAT_H_COMPLEX *gdualf);

extern void
LTFAT_H_NAME(gabtight_fac)(const LTFAT_H_COMPLEX *gf, const int L, 
			   const int a, const int M,
			   LTFAT_H_COMPLEX *gtightf);


void LTFAT_H_NAME(gabdual_long)(const LTFAT_H_COMPLEX *g,
				const int L, const int a,
				const int M, LTFAT_H_COMPLEX *gd);

void LTFAT_H_NAME(gabdualreal_long)(const LTFAT_H_REAL *g,
				    const int L, const int a,
				    const int M, LTFAT_H_REAL *gd);

void LTFAT_H_NAME(gabtight_long)(const LTFAT_H_COMPLEX *g,
				const int L, const int a,
				const int M, LTFAT_H_COMPLEX *gd);

void LTFAT_H_NAME(gabtightreal_long)(const LTFAT_H_REAL *g,
				    const int L, const int a,
				    const int M, LTFAT_H_REAL *gd);


/* --------- filter bank DGTs ---------------- */

extern void
LTFAT_H_NAME(dgt_fb)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
		     const int L, const int Lg,
		     const int W,  const int a, const int M, 
		     LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(dgt_fb_r)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		       const int L, const int gl,
		       const int W,  const int a, const int M, 
		       LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(dgtreal_fb)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		       const int L, const int gl,
		       const int W,  const int a, const int M, 
		       LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(idgt_fb)(const LTFAT_H_COMPLEX *cin, const LTFAT_H_COMPLEX *g,
		      const int L, const int gl,
		      const int W, const int a, const int M, 
		      LTFAT_H_COMPLEX *f);

/* -------- windows ------------------------------ */

extern void
LTFAT_H_NAME(pgauss)(const int L, const double w, const double c_t,
		     LTFAT_H_REAL *g);

extern void
LTFAT_H_NAME(pgauss_cmplx)(const int L, const double w, const double c_t, 
		     const double c_f, LTFAT_H_COMPLEX *g);

/* --------- other stuff -------- */
extern void
LTFAT_H_NAME(heapint)(const LTFAT_H_REAL *s,
		      const LTFAT_H_REAL *tgrad,
		      const LTFAT_H_REAL *fgrad,
		      const int a, const int M, const int L, const int W,
		      const LTFAT_H_REAL tol, LTFAT_H_REAL *phase);
		     
extern void
LTFAT_H_NAME(col2diag)(const LTFAT_H_COMPLEX *cin, const int L,
		       LTFAT_H_COMPLEX *cout);

extern void
LTFAT_H_NAME(col2diag_r)(const LTFAT_H_REAL *cin, const int L,
			 LTFAT_H_REAL *cout);

extern void
LTFAT_H_NAME(gabreassign)(const LTFAT_H_REAL *s, const LTFAT_H_REAL *tgrad,
			  const LTFAT_H_REAL *fgrad, const int L, const int W, 
			  const int a, const int M, LTFAT_H_REAL *sr);

extern void
LTFAT_H_NAME(fftshift_r)(const double *f, const int L, double *h);

extern void
LTFAT_H_NAME(ifftshift_r)(const LTFAT_H_REAL *f, const int L, LTFAT_H_REAL *h);

extern void
LTFAT_H_NAME(fir2iir_r)(const LTFAT_H_REAL *f, const int Lfir, const int Liir,
			LTFAT_H_REAL *h);
extern void
LTFAT_H_NAME(iir2fir_r)(const LTFAT_H_REAL *f, const int Liir,
			const int Lfir,	const int symm, LTFAT_H_REAL *h);

int
LTFAT_H_NAME(complexprod)(LTFAT_H_COMPLEX *c, const LTFAT_H_COMPLEX a,
			  const LTFAT_H_COMPLEX b);

/* ----- internal routines for calling BLAS and LAPACK ----- */

/* LAPACK overwrites the input argument. */
int
LTFAT_H_NAME(ltfat_posv)(const int N, const int NRHS,
			 LTFAT_H_COMPLEX *A, const int lda,
			 LTFAT_H_COMPLEX *B, const int ldb);

/* LAPACK overwrites the input argument. */
int
LTFAT_H_NAME(ltfat_gesvd)(const int M, const int N,
			  LTFAT_H_COMPLEX *A, const int lda,
			  LTFAT_H_REAL *S, LTFAT_H_COMPLEX *U, const int ldu,
			  LTFAT_H_COMPLEX *VT, const int ldvt);

void
LTFAT_H_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_TRANSPOSE TransB,
			 const int M, const int N, const int K,
			 const LTFAT_H_COMPLEX *alpha,
			 const LTFAT_H_COMPLEX *A, const int lda,
			 const LTFAT_H_COMPLEX *B, const int ldb,
			 const LTFAT_H_COMPLEX *beta,
			 LTFAT_H_COMPLEX *C, const int ldc);


/*  ----- experimental FFTW plan interface ------ */

typedef struct
{
  int a;
  int M;
  int L;
  int W;
  int c;
  int d;
  int h_a;
  LTFAT_H_FFTW(plan) p_before; 
  LTFAT_H_FFTW(plan) p_after;
  LTFAT_H_FFTW(plan) p_veryend;
  LTFAT_H_REAL *sbuf;
  const LTFAT_H_COMPLEX *f;
  LTFAT_H_COMPLEX *gf;
  LTFAT_H_COMPLEX *cout;
} LTFAT_H_NAME(ltfat_plan);

extern LTFAT_H_NAME(ltfat_plan)
LTFAT_H_NAME(plan_dgt_long)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
		       const int L, const int W, const int a,
		       const int M, LTFAT_H_COMPLEX *cout,
		       unsigned flags);

extern void 
LTFAT_H_NAME(ltfat_execute_plan)(const LTFAT_H_NAME(ltfat_plan) plan);

extern void
LTFAT_H_NAME(ltfat_destroy_plan)(LTFAT_H_NAME(ltfat_plan) plan);

extern void
LTFAT_H_NAME(dgt_walnut_plan)(LTFAT_H_NAME(ltfat_plan) plan,
			      const LTFAT_H_COMPLEX *f, 
			      const LTFAT_H_COMPLEX *gf);

