#include "dgt_long.h"
#include "dgt_multi.h"
#include "dgt_shear.h"

/*  --------- factorizations --------------- */

LTFAT_EXTERN void
LTFAT_H_NAME(wfac)(const LTFAT_H_COMPLEX *g, const int L, const int R,
		 const int a, const int M, LTFAT_H_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_H_NAME(wfac_r)(const LTFAT_H_REAL *g, const int L, const int R,
		   const int a, const int M, LTFAT_H_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_H_NAME(wfacreal)(const LTFAT_H_REAL *g, const int L, const int R,
		       const int a, const int M, LTFAT_H_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_H_NAME(iwfac)(const LTFAT_H_COMPLEX *gf, const int L, const int R,
		  const int a, const int M, LTFAT_H_COMPLEX *g);

LTFAT_EXTERN void
LTFAT_H_NAME(iwfac_r)(const LTFAT_H_COMPLEX *gf, const int L, const int R,
		    const int a, const int M, LTFAT_H_REAL *g);

LTFAT_EXTERN void
LTFAT_H_NAME(iwfacreal)(const LTFAT_H_COMPLEX *gf, const int L, const int R,
		    const int a, const int M, LTFAT_H_REAL *g);

/* --------- DGT by factorization ------------ */

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_fac)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *gf,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_H_COMPLEX *cout);


LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_long)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_fac_r)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf,
			const int L,
			const int W,  const int a,
			const int M, LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN
void LTFAT_H_NAME(dgtreal_fac)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf,
			       const int L,
			       const int W,  const int a,
			       const int M, LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_walnut_r)(const LTFAT_H_REAL *f, const LTFAT_H_COMPLEX *gf,
			 const int L, const int W,
			  const int a, const int M, LTFAT_H_COMPLEX *cout);


LTFAT_EXTERN void
LTFAT_H_NAME(idgt_fac)(const LTFAT_H_COMPLEX *c, const LTFAT_H_COMPLEX *gf,
		       const int L,
		       const int W,const int a, const int M,
		       LTFAT_H_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_H_NAME(idgtreal_fac)(const LTFAT_H_COMPLEX *cin, const LTFAT_H_COMPLEX *gf,
			   const int L, const int W,
			   const int a, const int M,
			   LTFAT_H_REAL *f);


/* --------- Wilson and WMDCT bases ---------*/
LTFAT_EXTERN void
LTFAT_H_NAME(dwilt_long)(const LTFAT_H_COMPLEX *f,
			     const LTFAT_H_COMPLEX *g,
			     const int L, const int W, const int M,
			     LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dwiltreal_long)(const LTFAT_H_REAL *f,
			     const LTFAT_H_REAL *g,
			     const int L, const int W, const int M,
			     LTFAT_H_REAL *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(wmdct_long)(const LTFAT_H_COMPLEX *f,
			     const LTFAT_H_COMPLEX *g,
			     const int L, const int W, const int M,
			     LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(wmdctreal_long)(const LTFAT_H_REAL *f,
			     const LTFAT_H_REAL *g,
			     const int L, const int W, const int M,
			     LTFAT_H_REAL *cout);

/* --------- dual windows etc. --------------- */

LTFAT_EXTERN void
LTFAT_H_NAME(gabdual_fac)(const LTFAT_H_COMPLEX *g, const int L, const int R,
			const int a, const int M, LTFAT_H_COMPLEX *gdualf);

LTFAT_EXTERN void
LTFAT_H_NAME(gabdualreal_fac)(const LTFAT_H_COMPLEX *g, const int L, const int R,
			const int a, const int M, LTFAT_H_COMPLEX *gdualf);

LTFAT_EXTERN void
LTFAT_H_NAME(gabtight_fac)(const LTFAT_H_COMPLEX *gf, const int L, const int R,
			   const int a, const int M,
			   LTFAT_H_COMPLEX *gtightf);

LTFAT_EXTERN void
LTFAT_H_NAME(gabtightreal_fac)(const LTFAT_H_COMPLEX *gf, const int L, const int R,
			   const int a, const int M,
			   LTFAT_H_COMPLEX *gtightf);


/* --------- filter bank DGTs ---------------- */

LTFAT_EXTERN void
LTFAT_H_NAME_COMPLEX(dgt_fb)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
		     const int L, const int Lg,
		     const int W,  const int a, const int M,
		     LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_fb)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		       const int L, const int gl,
		       const int W,  const int a, const int M,
		       LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_fb)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		       const int L, const int gl,
		       const int W,  const int a, const int M,
		       LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(idgt_fb)(const LTFAT_H_COMPLEX *cin, const LTFAT_H_COMPLEX *g,
		      const int L, const int gl,
		      const int W, const int a, const int M,
		      LTFAT_H_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_H_NAME(idgt_fb_r)(const LTFAT_H_COMPLEX *cin, const LTFAT_H_REAL *g,
		      const int L, const int gl,
		      const int W, const int a, const int M,
		      LTFAT_H_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_H_NAME(idgtreal_fb)(const LTFAT_H_COMPLEX *cin, const LTFAT_H_REAL *g,
			  const int L, const int gl, const int W,
			  const int a, const int M,
			  LTFAT_H_REAL *f);

/* ---------- OLA DGTs ------------- */
LTFAT_EXTERN void
LTFAT_H_NAME(dgt_ola)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
			  const int L, const int gl,
			  const int W, const int a, const int M, const int bl,
			  LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_ola)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
			  const int L, const int gl,
			  const int W, const int a, const int M, const int bl,
			  LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shearola)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
			 const int L, const int gl, const int W, const int a, const int M,
			 const int s0, const int s1, const int br, const int bl,
			 LTFAT_H_COMPLEX *cout);


/* --------- filterbank codes ------------*/
LTFAT_EXTERN void
LTFAT_H_NAME(ufilterbank_fft)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
                              const int L, const int gl,
			      const int W, const int a, const int M,
			      LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(convsub_fft)(const LTFAT_H_COMPLEXH *F, const LTFAT_H_COMPLEXH *G,
                          const size_t L, const size_t a, LTFAT_H_COMPLEXH *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(convsub_fft_plan)(const LTFAT_H_COMPLEXH *F, const LTFAT_H_COMPLEXH *G,
                               const size_t L, const size_t a, LTFAT_H_COMPLEXH *cout,
                               LTFAT_H_FFTW(plan) *p);

LTFAT_EXTERN void
LTFAT_H_NAME(convsub_fftbl)(const LTFAT_H_COMPLEXH *F, const size_t L,
                            const LTFAT_H_COMPLEXH *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_H_COMPLEXH *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(convsub_fftbl_plan)(const LTFAT_H_COMPLEXH *F, const size_t L,
                            const LTFAT_H_COMPLEXH *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_H_COMPLEXH *cout,
                            LTFAT_H_FFTW(plan) *p);

// Inverse
LTFAT_EXTERN void
LTFAT_H_NAME(upconv_fft)(const LTFAT_H_COMPLEXH *c, const size_t Lc,
                         const LTFAT_H_COMPLEXH *G, const size_t a,
                         LTFAT_H_COMPLEXH *Fout);

LTFAT_EXTERN void
LTFAT_H_NAME(upconv_fft_plan)(const LTFAT_H_COMPLEXH *c, const size_t Lc,
                              const LTFAT_H_COMPLEXH *G, const size_t a,
                              LTFAT_H_COMPLEXH *Fout, LTFAT_H_FFTW(plan) *p,
                              LTFAT_H_COMPLEXH *cbuffer);

LTFAT_EXTERN void
LTFAT_H_NAME(upconv_fftbl)(const LTFAT_H_COMPLEXH *c, const size_t Lc,
                            const LTFAT_H_COMPLEXH *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_H_COMPLEXH *Fout);

LTFAT_EXTERN void
LTFAT_H_NAME(upconv_fftbl_plan)(const LTFAT_H_COMPLEXH *c, const size_t Lc,
                            const LTFAT_H_COMPLEXH *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_H_COMPLEXH *Fout,
                            LTFAT_H_FFTW(plan) *p,LTFAT_H_COMPLEXH *cbuffer);



/* -------- windows ------------------------------ */

LTFAT_EXTERN void
LTFAT_H_NAME(pgauss)(const int L, const double w, const double c_t,
		     LTFAT_H_REAL *g);

LTFAT_EXTERN void
LTFAT_H_NAME(pgauss_cmplx)(const int L, const double w, const double c_t,
		     const double c_f, LTFAT_H_COMPLEX *g);


/* --------- pfilt and filterbanks ------------- */
LTFAT_EXTERN void
LTFAT_H_NAME(pfilt_fir_rr)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
			   const int L, const int gl,
			   const int W, const int a,
			   LTFAT_H_REAL *cout);

/* --------- other stuff -------- */
LTFAT_EXTERN void
LTFAT_H_NAME(heapint)(const LTFAT_H_REAL *s,
		      const LTFAT_H_REAL *tgrad,
		      const LTFAT_H_REAL *fgrad,
		      const int a, const int M, const int L, const int W,
		      const LTFAT_H_REAL tol, LTFAT_H_REAL *phase);

LTFAT_EXTERN void
LTFAT_H_NAME(gabreassign)(const LTFAT_H_REAL *s, const LTFAT_H_REAL *tgrad,
			  const LTFAT_H_REAL *fgrad, const int L, const int W,
			  const int a, const int M, LTFAT_H_REAL *sr);

LTFAT_EXTERN void
LTFAT_H_NAME(fftshift_r)(const LTFAT_H_REAL *f, const int L, LTFAT_H_REAL *h);

LTFAT_EXTERN void
LTFAT_H_NAME(ifftshift_r)(const LTFAT_H_REAL *f, const int L, LTFAT_H_REAL *h);

LTFAT_EXTERN void
LTFAT_H_NAME(fir2long_r)(const LTFAT_H_REAL *f, const int Lfir, const int Llong,
			LTFAT_H_REAL *h);

LTFAT_EXTERN void
LTFAT_H_NAME(fir2long_c)(const LTFAT_H_COMPLEX *f,
			 const int Lfir, const int Llong,
			 LTFAT_H_COMPLEX *h);

LTFAT_EXTERN void
LTFAT_H_NAME(long2fir_r)(const LTFAT_H_REAL *f, const int Llong,
			const int Lfir, LTFAT_H_REAL *h);

LTFAT_EXTERN void
LTFAT_H_NAME(long2fir_c)(const LTFAT_H_COMPLEX *f, const int Llong,
			const int Lfir, LTFAT_H_COMPLEX *h);

int
LTFAT_H_NAME(complexprod)(LTFAT_H_COMPLEX *c, const LTFAT_H_COMPLEX a,
			  const LTFAT_H_COMPLEX b);

/* ----- internal routines for calling BLAS and LAPACK ----- */

/*
// LAPACK overwrites the input argument.
int
LTFAT_H_NAME(ltfat_posv)(const int N, const int NRHS,
			 LTFAT_H_COMPLEX *A, const int lda,
			 LTFAT_H_COMPLEX *B, const int ldb);

// LAPACK overwrites the input argument.
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
*/


// LAPACK overwrites the input argument.
int
LTFAT_H_NAME(ltfat_posv)(const ptrdiff_t N, const ptrdiff_t NRHS,
			 LTFAT_H_COMPLEX *A, const ptrdiff_t lda,
			 LTFAT_H_COMPLEX *B, const ptrdiff_t ldb);

// LAPACK overwrites the input argument.
int
LTFAT_H_NAME(ltfat_gesvd)(const ptrdiff_t M, const ptrdiff_t N,
			  LTFAT_H_COMPLEX *A, const ptrdiff_t lda,
			  LTFAT_H_REAL *S, LTFAT_H_COMPLEX *U, const ptrdiff_t ldu,
			  LTFAT_H_COMPLEX *VT, const ptrdiff_t ldvt);

void
LTFAT_H_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_TRANSPOSE TransB,
			 const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
			 const LTFAT_H_COMPLEX *alpha,
			 const LTFAT_H_COMPLEX *A, const ptrdiff_t lda,
			 const LTFAT_H_COMPLEX *B, const ptrdiff_t ldb,
			 const LTFAT_H_COMPLEX *beta,
			 LTFAT_H_COMPLEX *C, const ptrdiff_t ldc);

/*   --- dgtreal_long class definition  --- */
typedef struct
{
  int a;
  int M;
  int L;
  int W;
  int c;
  int h_a;
  LTFAT_H_FFTW(plan) p_before;
  LTFAT_H_FFTW(plan) p_after;
  LTFAT_H_FFTW(plan) p_veryend;
  LTFAT_H_REAL *sbuf;
  LTFAT_H_COMPLEX *cbuf;
  const LTFAT_H_REAL *f;
  LTFAT_H_COMPLEX *gf;
  LTFAT_H_REAL *cwork;
  LTFAT_H_COMPLEX *cout;
  LTFAT_H_REAL *ff, *cf;
} LTFAT_H_NAME(dgtreal_long_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgtreal_long_plan)
LTFAT_H_NAME(dgtreal_long_init)(const LTFAT_H_REAL *f, const LTFAT_H_REAL *g,
		       const int L, const int W, const int a,
		       const int M, LTFAT_H_COMPLEX *cout,
		       unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_long_execute)(const LTFAT_H_NAME(dgtreal_long_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_long_done)(LTFAT_H_NAME(dgtreal_long_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_walnut_plan)(LTFAT_H_NAME(dgtreal_long_plan) plan);




/*   --- dgt_fb class definition  --- */
typedef struct
{
  int a;
  int M;
  int gl;

  LTFAT_H_FFTW(plan) p_small;
  LTFAT_H_REAL *sbuf;
  LTFAT_H_REAL *fw;
  LTFAT_H_COMPLEX *gw;
  LTFAT_H_COMPLEX *cout;
} LTFAT_H_NAME(dgt_fb_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgt_fb_plan)
LTFAT_H_NAME(dgt_fb_init)(const LTFAT_H_COMPLEX *g,
   const int gl, const int a, const int M,
   unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_fb_execute)(const LTFAT_H_NAME(dgt_fb_plan) plan,
			     const LTFAT_H_COMPLEX *f, const int L, const int W,
			     LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_fb_done)(LTFAT_H_NAME(dgt_fb_plan) plan);



/*   --- dgtreal_fb class definition  --- */

typedef struct
{
  int a;
  int M;
  int gl;

  LTFAT_H_FFTW(plan) p_small;
  LTFAT_H_REAL    *sbuf;
  LTFAT_H_COMPLEX *cbuf;
  LTFAT_H_REAL *fw;
  LTFAT_H_REAL *gw;
  LTFAT_H_COMPLEX *cout;
} LTFAT_H_NAME(dgtreal_fb_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgtreal_fb_plan)
LTFAT_H_NAME(dgtreal_fb_init)(const LTFAT_H_REAL *g,
   const int gl, const int a, const int M,
   unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_fb_execute)(const LTFAT_H_NAME(dgtreal_fb_plan) plan,
			     const LTFAT_H_REAL *f, const int L, const int W,
			     LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_fb_done)(LTFAT_H_NAME(dgtreal_fb_plan) plan);


/*   --- dgt_ola class definition  --- */
typedef struct
{
   LTFAT_H_NAME(dgt_long_plan) plan;
   int bl;
   int gl;
   int W;
   LTFAT_H_COMPLEX *buf;
   LTFAT_H_COMPLEX *gext;
   LTFAT_H_COMPLEX *cbuf;

} LTFAT_H_NAME(dgt_ola_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgt_ola_plan)
LTFAT_H_NAME(dgt_ola_init)(const LTFAT_H_COMPLEX *g, const int gl,
			   const int W, const int a, const int M, const int bl,
			   unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_ola_execute)(const LTFAT_H_NAME(dgt_ola_plan) plan,
			    const LTFAT_H_COMPLEX *f, const int L,
			    LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_ola_done)(LTFAT_H_NAME(dgt_ola_plan) plan);


LTFAT_EXTERN void
LTFAT_H_NAME(dgt_walnut_plan)(LTFAT_H_NAME(dgt_long_plan) plan);


/*   --- dgtreal_ola class definition  --- */
typedef struct
{
   LTFAT_H_NAME(dgtreal_long_plan) plan;
   int bl;
   int gl;
   int W;
   LTFAT_H_REAL *buf;
   LTFAT_H_REAL *gext;
   LTFAT_H_COMPLEX *cbuf;

} LTFAT_H_NAME(dgtreal_ola_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgtreal_ola_plan)
LTFAT_H_NAME(dgtreal_ola_init)(const LTFAT_H_REAL *g, const int gl,
			       const int W, const int a, const int M, const int bl,
			       unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_ola_execute)(const LTFAT_H_NAME(dgtreal_ola_plan) plan,
				  const LTFAT_H_REAL *f, const int L,
				  LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_ola_done)(LTFAT_H_NAME(dgtreal_ola_plan) plan);


LTFAT_EXTERN void
LTFAT_H_NAME(dgtreal_walnut_plan)(LTFAT_H_NAME(dgtreal_long_plan) plan);

/* -----  dgt_shearola class definition ------ */

typedef struct
{
   LTFAT_H_NAME(dgt_shear_plan) plan;
   int bl;
   int gl;
   int W;
   LTFAT_H_COMPLEX *buf;
   LTFAT_H_COMPLEX *gext;
   LTFAT_H_COMPLEX *cbuf;

} LTFAT_H_NAME(dgt_shearola_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgt_shearola_plan)
LTFAT_H_NAME(dgt_shearola_init)(const LTFAT_H_COMPLEX *g, const int gl,
			   const int W, const int a, const int M,
			   const int s0, const int s1, const int br,
			   const int bl,
			   unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shearola_execute)(const LTFAT_H_NAME(dgt_shearola_plan) plan,
			    const LTFAT_H_COMPLEX *f, const int L,
			    LTFAT_H_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shearola_done)(LTFAT_H_NAME(dgt_shearola_plan) plan);

