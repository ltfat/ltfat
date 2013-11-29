#include "dgt_long.h"
#include "dgt_multi.h"
#include "dgt_shear.h"

/*  --------- factorizations --------------- */

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(wfac)(const LTFAT_COMPLEX *g, const int L, const int R,
		 const int a, const int M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME_REAL(wfac)(const LTFAT_REAL *g, const int L, const int R,
		   const int a, const int M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME(wfacreal)(const LTFAT_REAL *g, const int L, const int R,
		       const int a, const int M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(iwfac)(const LTFAT_COMPLEX *gf, const int L, const int R,
		  const int a, const int M, LTFAT_COMPLEX *g);

LTFAT_EXTERN void
LTFAT_NAME_REAL(iwfac)(const LTFAT_COMPLEX *gf, const int L, const int R,
		    const int a, const int M, LTFAT_REAL *g);

LTFAT_EXTERN void
LTFAT_NAME(iwfacreal)(const LTFAT_COMPLEX *gf, const int L, const int R,
		    const int a, const int M, LTFAT_REAL *g);

/* --------- DGT by factorization ------------ */

LTFAT_EXTERN void
LTFAT_NAME(dgt_fac)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *gf,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_COMPLEX *cout);


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
		      const int L, const int W,  const int a,
		      const int M, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fac_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
			const int L,
			const int W,  const int a,
			const int M, LTFAT_COMPLEX *cout);

LTFAT_EXTERN
void LTFAT_NAME(dgtreal_fac)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
			       const int L,
			       const int W,  const int a,
			       const int M, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_walnut_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
			 const int L, const int W,
			  const int a, const int M, LTFAT_COMPLEX *cout);


LTFAT_EXTERN void
LTFAT_NAME(idgt_fac)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *gf,
		       const int L,
		       const int W,const int a, const int M,
		       LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf,
			   const int L, const int W,
			   const int a, const int M,
			   LTFAT_REAL *f);

LTFAT_EXTERN void
LTFAT_NAME(idgt_long)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
			  const int L, const int W,
			  const int a, const int M,
			  LTFAT_COMPLEX *f);


LTFAT_EXTERN void
LTFAT_NAME(idgtreal_long)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
			      const int L, const int W,
			      const int a, const int M,
			      LTFAT_REAL *f);


/* --------- dual windows etc. --------------- */

LTFAT_EXTERN void
LTFAT_NAME(gabdual_fac)(const LTFAT_COMPLEX *g, const int L, const int R,
			const int a, const int M, LTFAT_COMPLEX *gdualf);

LTFAT_EXTERN void
LTFAT_NAME(gabdualreal_fac)(const LTFAT_COMPLEX *g, const int L, const int R,
			const int a, const int M, LTFAT_COMPLEX *gdualf);

LTFAT_EXTERN void
LTFAT_NAME(gabtight_fac)(const LTFAT_COMPLEX *gf, const int L, const int R,
			   const int a, const int M,
			   LTFAT_COMPLEX *gtightf);

LTFAT_EXTERN void
LTFAT_NAME(gabtightreal_fac)(const LTFAT_COMPLEX *gf, const int L, const int R,
			   const int a, const int M,
			   LTFAT_COMPLEX *gtightf);


/* --------- filter bank DGTs ---------------- */

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dgt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
		     const int L, const int gl,
		     const int W,  const int a, const int M,
		     LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
		       const int L, const int gl,
		       const int W,  const int a, const int M,
		       LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
		       const int L, const int gl,
		       const int W,  const int a, const int M,
		       LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
		      const int L, const int gl,
		      const int W, const int a, const int M,
		      LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb_r)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
		      const int L, const int gl,
		      const int W, const int a, const int M,
		      LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
			  const int L, const int gl, const int W,
			  const int a, const int M,
			  LTFAT_REAL *f);

/* ---------- OLA DGTs ------------- */
LTFAT_EXTERN void
LTFAT_NAME(dgt_ola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			  const int L, const int gl,
			  const int W, const int a, const int M, const int bl,
			  LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			  const int L, const int gl,
			  const int W, const int a, const int M, const int bl,
			  LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			 const int L, const int gl, const int W, const int a, const int M,
			 const int s0, const int s1, const int br, const int bl,
			 LTFAT_COMPLEX *cout);


/* --------- filterbank codes ------------*/
LTFAT_EXTERN void
LTFAT_NAME(ufilterbank_fft)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                              const int L, const int gl,
			      const int W, const int a, const int M,
			      LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                          const size_t L, const size_t a, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_plan)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                               const size_t L, const size_t a, LTFAT_COMPLEX *cout,
                               LTFAT_FFTW(plan) *p);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl)(const LTFAT_COMPLEX *F, const size_t L,
                            const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_plan)(const LTFAT_COMPLEX *F, const size_t L,
                            const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_COMPLEX *cout,
                            LTFAT_FFTW(plan) *p);

// Inverse
LTFAT_EXTERN void
LTFAT_NAME(upconv_fft)(const LTFAT_COMPLEX *c, const size_t Lc,
                         const LTFAT_COMPLEX *G, const size_t a,
                         LTFAT_COMPLEX *Fout);

LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_plan)(const LTFAT_COMPLEX *c, const size_t Lc,
                              const LTFAT_COMPLEX *G, const size_t a,
                              LTFAT_COMPLEX *Fout, LTFAT_FFTW(plan) *p,
                              LTFAT_COMPLEX *cbuffer);

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl)(const LTFAT_COMPLEX *c, const size_t Lc,
                            const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_COMPLEX *Fout);

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_plan)(const LTFAT_COMPLEX *c, const size_t Lc,
                            const LTFAT_COMPLEX *G, const size_t Lg, const ptrdiff_t foff,
                            const double a, const double realonly, LTFAT_COMPLEX *Fout,
                            LTFAT_FFTW(plan) *p,LTFAT_COMPLEX *cbuffer);



/* -------- windows ------------------------------ */

LTFAT_EXTERN void
LTFAT_NAME(pgauss)(const int L, const double w, const double c_t,
		     LTFAT_REAL *g);

LTFAT_EXTERN void
LTFAT_NAME(pgauss_cmplx)(const int L, const double w, const double c_t,
		     const double c_f, LTFAT_COMPLEX *g);


/* --------- pfilt and filterbanks ------------- */
LTFAT_EXTERN void
LTFAT_NAME(pfilt_fir_rr)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			   const int L, const int gl,
			   const int W, const int a,
			   LTFAT_REAL *cout);

/* --------- other stuff -------- */
LTFAT_EXTERN void
LTFAT_NAME(heapint)(const LTFAT_REAL *s,
		      const LTFAT_REAL *tgrad,
		      const LTFAT_REAL *fgrad,
		      const int a, const int M, const int L, const int W,
		      const LTFAT_REAL tol, LTFAT_REAL *phase);

LTFAT_EXTERN void
LTFAT_NAME(gabreassign)(const LTFAT_REAL *s, const LTFAT_REAL *tgrad,
			  const LTFAT_REAL *fgrad, const int L, const int W,
			  const int a, const int M, LTFAT_REAL *sr);

LTFAT_EXTERN void
LTFAT_NAME(fftshift_r)(const LTFAT_REAL *f, const int L, LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(ifftshift_r)(const LTFAT_REAL *f, const int L, LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(fir2long_r)(const LTFAT_REAL *f, const int Lfir, const int Llong,
			LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(fir2long_c)(const LTFAT_COMPLEX *f,
			 const int Lfir, const int Llong,
			 LTFAT_COMPLEX *h);

LTFAT_EXTERN void
LTFAT_NAME(long2fir_r)(const LTFAT_REAL *f, const int Llong,
			const int Lfir, LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(long2fir_c)(const LTFAT_COMPLEX *f, const int Llong,
			const int Lfir, LTFAT_COMPLEX *h);

int
LTFAT_NAME(complexprod)(LTFAT_COMPLEX *c, const LTFAT_COMPLEX a,
			  const LTFAT_COMPLEX b);

/* ----- internal routines for calling BLAS and LAPACK ----- */

/*
// LAPACK overwrites the input argument.
int
LTFAT_NAME(ltfat_posv)(const int N, const int NRHS,
			 LTFAT_COMPLEX *A, const int lda,
			 LTFAT_COMPLEX *B, const int ldb);

// LAPACK overwrites the input argument.
int
LTFAT_NAME(ltfat_gesvd)(const int M, const int N,
			  LTFAT_COMPLEX *A, const int lda,
			  LTFAT_REAL *S, LTFAT_COMPLEX *U, const int ldu,
			  LTFAT_COMPLEX *VT, const int ldvt);

void
LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_TRANSPOSE TransB,
			 const int M, const int N, const int K,
			 const LTFAT_COMPLEX *alpha,
			 const LTFAT_COMPLEX *A, const int lda,
			 const LTFAT_COMPLEX *B, const int ldb,
			 const LTFAT_COMPLEX *beta,
			 LTFAT_COMPLEX *C, const int ldc);
*/


// LAPACK overwrites the input argument.
int
LTFAT_NAME(ltfat_posv)(const ptrdiff_t N, const ptrdiff_t NRHS,
			 LTFAT_COMPLEX *A, const ptrdiff_t lda,
			 LTFAT_COMPLEX *B, const ptrdiff_t ldb);

// LAPACK overwrites the input argument.
int
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

/*   --- dgtreal_long class definition  --- */
typedef struct
{
  int a;
  int M;
  int L;
  int W;
  int c;
  int h_a;
  LTFAT_FFTW(plan) p_before;
  LTFAT_FFTW(plan) p_after;
  LTFAT_FFTW(plan) p_veryend;
  LTFAT_REAL *sbuf;
  LTFAT_COMPLEX *cbuf;
  const LTFAT_REAL *f;
  LTFAT_COMPLEX *gf;
  LTFAT_REAL *cwork;
  LTFAT_COMPLEX *cout;
  LTFAT_REAL *ff, *cf;
} LTFAT_NAME(dgtreal_long_plan);


LTFAT_EXTERN LTFAT_NAME(dgtreal_long_plan)
LTFAT_NAME(dgtreal_long_init)(const LTFAT_REAL *f, const LTFAT_REAL *g,
		       const int L, const int W, const int a,
		       const int M, LTFAT_COMPLEX *cout,
		       unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long_execute)(const LTFAT_NAME(dgtreal_long_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan) plan);




/*   --- dgt_fb class definition  --- */
typedef struct
{
  int a;
  int M;
  int gl;

  LTFAT_FFTW(plan) p_small;
  LTFAT_REAL *sbuf;
  LTFAT_REAL *fw;
  LTFAT_COMPLEX *gw;
  LTFAT_COMPLEX *cout;
} LTFAT_NAME(dgt_fb_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_fb_plan)
LTFAT_NAME(dgt_fb_init)(const LTFAT_COMPLEX *g,
   const int gl, const int a, const int M,
   unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_execute)(const LTFAT_NAME(dgt_fb_plan) plan,
			     const LTFAT_COMPLEX *f, const int L, const int W,
			     LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_done)(LTFAT_NAME(dgt_fb_plan) plan);



/*   --- dgtreal_fb class definition  --- */

typedef struct
{
  int a;
  int M;
  int gl;

  LTFAT_FFTW(plan) p_small;
  LTFAT_REAL    *sbuf;
  LTFAT_COMPLEX *cbuf;
  LTFAT_REAL *fw;
  LTFAT_REAL *gw;
  LTFAT_COMPLEX *cout;
} LTFAT_NAME(dgtreal_fb_plan);


LTFAT_EXTERN LTFAT_NAME(dgtreal_fb_plan)
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL *g,
   const int gl, const int a, const int M,
   unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_execute)(const LTFAT_NAME(dgtreal_fb_plan) plan,
			     const LTFAT_REAL *f, const int L, const int W,
			     LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan) plan);


/*   --- dgt_ola class definition  --- */
typedef struct
{
   LTFAT_NAME(dgt_long_plan) plan;
   int bl;
   int gl;
   int W;
   LTFAT_COMPLEX *buf;
   LTFAT_COMPLEX *gext;
   LTFAT_COMPLEX *cbuf;

} LTFAT_NAME(dgt_ola_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_ola_plan)
LTFAT_NAME(dgt_ola_init)(const LTFAT_COMPLEX *g, const int gl,
			   const int W, const int a, const int M, const int bl,
			   unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_execute)(const LTFAT_NAME(dgt_ola_plan) plan,
			    const LTFAT_COMPLEX *f, const int L,
			    LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_done)(LTFAT_NAME(dgt_ola_plan) plan);


LTFAT_EXTERN void
LTFAT_NAME(dgt_walnut_plan)(LTFAT_NAME(dgt_long_plan) plan);


/*   --- dgtreal_ola class definition  --- */
typedef struct
{
   LTFAT_NAME(dgtreal_long_plan) plan;
   int bl;
   int gl;
   int W;
   LTFAT_REAL *buf;
   LTFAT_REAL *gext;
   LTFAT_COMPLEX *cbuf;

} LTFAT_NAME(dgtreal_ola_plan);


LTFAT_EXTERN LTFAT_NAME(dgtreal_ola_plan)
LTFAT_NAME(dgtreal_ola_init)(const LTFAT_REAL *g, const int gl,
			       const int W, const int a, const int M, const int bl,
			       unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_execute)(const LTFAT_NAME(dgtreal_ola_plan) plan,
				  const LTFAT_REAL *f, const int L,
				  LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_done)(LTFAT_NAME(dgtreal_ola_plan) plan);


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan) plan);

/* -----  dgt_shearola class definition ------ */

typedef struct
{
   LTFAT_NAME(dgt_shear_plan) plan;
   int bl;
   int gl;
   int W;
   LTFAT_COMPLEX *buf;
   LTFAT_COMPLEX *gext;
   LTFAT_COMPLEX *cbuf;

} LTFAT_NAME(dgt_shearola_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_shearola_plan)
LTFAT_NAME(dgt_shearola_init)(const LTFAT_COMPLEX *g, const int gl,
			   const int W, const int a, const int M,
			   const int s0, const int s1, const int br,
			   const int bl,
			   unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_execute)(const LTFAT_NAME(dgt_shearola_plan) plan,
			    const LTFAT_COMPLEX *f, const int L,
			    LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_done)(LTFAT_NAME(dgt_shearola_plan) plan);

