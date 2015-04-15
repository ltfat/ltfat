#include "dgt_long.h"
#include "dgt_multi.h"
#include "dgt_shear.h"

/*  --------- factorizations --------------- */

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(wfac)(const LTFAT_COMPLEX *g, const ltfatInt L, const ltfatInt R,
                         const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME_REAL(wfac)(const LTFAT_REAL *g, const ltfatInt L, const ltfatInt R,
                      const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME(wfacreal)(const LTFAT_REAL *g, const ltfatInt L, const ltfatInt R,
                     const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(iwfac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                          const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *g);

LTFAT_EXTERN void
LTFAT_NAME_REAL(iwfac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                       const ltfatInt a, const ltfatInt M, LTFAT_REAL *g);

LTFAT_EXTERN void
LTFAT_NAME(iwfacreal)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                      const ltfatInt a, const ltfatInt M, LTFAT_REAL *g);

/* --------- DGT by factorization ------------ */

LTFAT_EXTERN void
LTFAT_NAME(dgt_fac)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *gf,
                    const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                    const ltfatInt M, LTFAT_COMPLEX *cout);


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                         const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                         const ltfatInt M, const dgt_phasetype ptype,
                         LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fac_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
                      const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                      const ltfatInt M, const dgt_phasetype ptype,
                      LTFAT_COMPLEX *cout);

LTFAT_EXTERN
void LTFAT_NAME(dgtreal_fac)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
                             const ltfatInt L,
                             const ltfatInt W,  const ltfatInt a,
                             const ltfatInt M, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_walnut_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
                         const ltfatInt L, const ltfatInt W,
                         const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *cout);


LTFAT_EXTERN void
LTFAT_NAME(idgt_fac)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *gf,
                     const ltfatInt L,
                     const ltfatInt W, const ltfatInt a, const ltfatInt M,
                     const dgt_phasetype ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fac)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *gf,
                         const ltfatInt L, const ltfatInt W,
                         const ltfatInt a, const ltfatInt M,
                         const dgt_phasetype ptype, LTFAT_REAL *f);

LTFAT_EXTERN void
LTFAT_NAME(idgt_long)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
                      const ltfatInt L, const ltfatInt W,
                      const ltfatInt a, const ltfatInt M,
                      const dgt_phasetype ptype, LTFAT_COMPLEX *f);


LTFAT_EXTERN void
LTFAT_NAME(idgtreal_long)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M,
                          const dgt_phasetype ptype, LTFAT_REAL *f);


/* --------- dual windows etc. --------------- */

LTFAT_EXTERN void
LTFAT_NAME(gabdual_fac)(const LTFAT_COMPLEX *g, const ltfatInt L, const ltfatInt R,
                        const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *gdualf);

LTFAT_EXTERN void
LTFAT_NAME(gabdualreal_fac)(const LTFAT_COMPLEX *g, const ltfatInt L, const ltfatInt R,
                            const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *gdualf);

LTFAT_EXTERN void
LTFAT_NAME(gabtight_fac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                         const ltfatInt a, const ltfatInt M,
                         LTFAT_COMPLEX *gtightf);

LTFAT_EXTERN void
LTFAT_NAME(gabtightreal_fac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                             const ltfatInt a, const ltfatInt M,
                             LTFAT_COMPLEX *gtightf);


/* --------- filter bank DGTs ---------------- */

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dgt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                           const ltfatInt L, const ltfatInt gl,
                           const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                           const dgt_phasetype ptype, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                   const ltfatInt L, const ltfatInt gl,
                   const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                   const dgt_phasetype ptype, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                       const ltfatInt L, const ltfatInt gl,
                       const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                       const dgt_phasetype ptype, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
                    const dgt_phasetype ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgt_fb_r)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                      const ltfatInt L, const ltfatInt gl,
                      const ltfatInt W, const ltfatInt a, const ltfatInt M,
                      LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                        const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, LTFAT_REAL *f);

/* ---------- OLA DGTs ------------- */
LTFAT_EXTERN void
LTFAT_NAME(dgt_ola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a,
                    const ltfatInt M, const ltfatInt bl,
                    const dgt_phasetype ptype,
                    LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                        const ltfatInt L, const ltfatInt gl,
                        const ltfatInt W, const ltfatInt a, const ltfatInt M,
                        const ltfatInt bl, const dgt_phasetype ptype,
                        LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                         const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                         const ltfatInt a, const ltfatInt M,
                         const ltfatInt s0, const ltfatInt s1, const ltfatInt br,
                         const ltfatInt bl, LTFAT_COMPLEX *cout);

/* --------- FFT ------------------*/
LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(fftreal_init)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
                         LTFAT_COMPLEX *cout, unsigned flag);

LTFAT_EXTERN void
LTFAT_NAME(fftreal_execute)(LTFAT_FFTW(plan) p, LTFAT_REAL *f,
                            LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(fftreal)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
                    LTFAT_COMPLEX *cout);

LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(ifftreal_init)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
                          LTFAT_REAL *f, unsigned flag);

LTFAT_EXTERN void
LTFAT_NAME(ifftreal_execute)(LTFAT_FFTW(plan), LTFAT_COMPLEX *c,
                             const ltfatInt L, const ltfatInt W,
                             LTFAT_REAL *f);

LTFAT_EXTERN void
LTFAT_NAME(ifftreal)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
                     LTFAT_REAL *f);


/* --------- filterbank codes ------------*/

typedef struct LTFAT_NAME(convsub_fft_plan_struct) *LTFAT_NAME(convsub_fft_plan);
typedef struct LTFAT_NAME(convsub_fftbl_plan_struct) *LTFAT_NAME(convsub_fftbl_plan);

typedef struct LTFAT_NAME(upconv_fft_plan_struct) *LTFAT_NAME(upconv_fft_plan);
typedef struct LTFAT_NAME(upconv_fftbl_plan_struct) *LTFAT_NAME(upconv_fftbl_plan);

LTFAT_EXTERN void
LTFAT_NAME(ufilterbank_fft)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                            const ltfatInt L, const ltfatInt gl,
                            const ltfatInt W, const ltfatInt a, const ltfatInt M,
                            LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                           const ltfatInt L, const ltfatInt W,
                           const ltfatInt a[], const ltfatInt M,
                           LTFAT_COMPLEX *cout[]);

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fft_execute)(LTFAT_NAME(convsub_fft_plan) p[],
                                   const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                                   const ltfatInt M, LTFAT_COMPLEX *cout[]);


LTFAT_EXTERN LTFAT_NAME(convsub_fft_plan)
LTFAT_NAME(convsub_fft_init)(const ltfatInt L, const ltfatInt W,
                             const ltfatInt a, const LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_done)(LTFAT_NAME(convsub_fft_plan) p);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft_execute)(const LTFAT_NAME(convsub_fft_plan) p,
                                const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                                LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fft)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G,
                        const ltfatInt L, const ltfatInt W, const ltfatInt a,
                        LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fftbl)(const LTFAT_COMPLEX *F, const LTFAT_COMPLEX *G[],
                             const ltfatInt L, const ltfatInt Gl[],
                             const ltfatInt W, const double a[], const ltfatInt M,
                             const ltfatInt foff[], const int realonly[],
                             LTFAT_COMPLEX *cout[]);

LTFAT_EXTERN void
LTFAT_NAME(filterbank_fftbl_execute)(LTFAT_NAME(convsub_fftbl_plan) p[],
                                     const LTFAT_COMPLEX *F,
                                     const LTFAT_COMPLEX *G[],
                                     const ltfatInt M, const ltfatInt foff[],
                                     const int realonly[], LTFAT_COMPLEX *cout[]);

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft_execute)(LTFAT_NAME(upconv_fft_plan) p[],
                                    const LTFAT_COMPLEX *cin[],
                                    const LTFAT_COMPLEX *G[],
                                    const ltfatInt M,
                                    LTFAT_COMPLEX *F );


LTFAT_EXTERN LTFAT_NAME(convsub_fftbl_plan)
LTFAT_NAME(convsub_fftbl_init)( const ltfatInt L, const ltfatInt Gl,
                                const ltfatInt W, const double a,
                                const LTFAT_COMPLEX *cout);

LTFAT_EXTERN LTFAT_NAME(convsub_fftbl_plan)
LTFAT_NAME(convsub_fftbl_init_no_ifft_plan)( const ltfatInt L, const ltfatInt Gl,
                                             const ltfatInt W, const double a);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_done)( LTFAT_NAME(convsub_fftbl_plan) p);


LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl_execute)(const LTFAT_NAME(convsub_fftbl_plan) p,
                                  const LTFAT_COMPLEX *F,
                                  const LTFAT_COMPLEX *G,
                                  const ltfatInt foff,
                                  const int realonly,
                                  LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(convsub_fftbl)(const LTFAT_COMPLEX *F,  const LTFAT_COMPLEX *G,
                          const ltfatInt L, const ltfatInt Gl, const ltfatInt W,
                          const double a, const ltfatInt foff,
                          const int realonly, LTFAT_COMPLEX *cout);



// Inverse
LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                            const ltfatInt L, const ltfatInt W,
                            const ltfatInt a[], const ltfatInt M,
                            LTFAT_COMPLEX *F);

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fft_execute)(LTFAT_NAME(upconv_fft_plan) p[],
                                    const LTFAT_COMPLEX *cin[],
                                    const LTFAT_COMPLEX *G[],
                                    const ltfatInt M,
                                    LTFAT_COMPLEX *F );


LTFAT_EXTERN void
LTFAT_NAME(upconv_fft)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                       const ltfatInt L, const ltfatInt W, const ltfatInt a,
                       LTFAT_COMPLEX *F);

LTFAT_EXTERN LTFAT_NAME(upconv_fft_plan)
LTFAT_NAME(upconv_fft_init)(const ltfatInt L, const ltfatInt W, const ltfatInt a);


LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_execute)(LTFAT_NAME(upconv_fft_plan) p,
                               const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                               LTFAT_COMPLEX *F);

LTFAT_EXTERN void
LTFAT_NAME(upconv_fft_done)(LTFAT_NAME(upconv_fft_plan) p);



LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fftbl)(const LTFAT_COMPLEX *cin[], const LTFAT_COMPLEX *G[],
                              const ltfatInt L, const ltfatInt Gl[],
                              const ltfatInt W, const double a[], const ltfatInt M,
                              const ptrdiff_t foff[], const int realonly[],
                              LTFAT_COMPLEX *F);

LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_fftbl_execute)(LTFAT_NAME(upconv_fftbl_plan) p[],
                                      const LTFAT_COMPLEX *cin[],
                                      const LTFAT_COMPLEX *G[],
                                      const ltfatInt M, const ltfatInt foff[],
                                      const int realonly[],
                                      LTFAT_COMPLEX *F);



LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *G,
                         const ltfatInt L, const ltfatInt Gl, const ltfatInt W,
                         const double a, const ptrdiff_t foff,
                         const int realonly, LTFAT_COMPLEX *F);

LTFAT_EXTERN LTFAT_NAME(upconv_fftbl_plan)
LTFAT_NAME(upconv_fftbl_init)( const ltfatInt L, const ltfatInt Gl,
                               const ltfatInt W, const double a);

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_done)(LTFAT_NAME(upconv_fftbl_plan) p);

LTFAT_EXTERN void
LTFAT_NAME(upconv_fftbl_execute)(const LTFAT_NAME(upconv_fftbl_plan) p,
                                 const LTFAT_COMPLEX *cin,
                                 const LTFAT_COMPLEX *G,
                                 const ptrdiff_t foff, const int realonly,
                                 LTFAT_COMPLEX *F);



/* -------- windows ------------------------------ */

LTFAT_EXTERN void
LTFAT_NAME(pgauss)(const ltfatInt L, const double w, const double c_t,
                   LTFAT_REAL *g);

LTFAT_EXTERN void
LTFAT_NAME(pgauss_cmplx)(const ltfatInt L, const double w, const double c_t,
                         const double c_f, LTFAT_COMPLEX *g);


/* --------- pfilt and filterbanks ------------- */
LTFAT_EXTERN void
LTFAT_NAME(pfilt_fir_rr)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                         const ltfatInt L, const ltfatInt gl,
                         const ltfatInt W, const ltfatInt a,
                         LTFAT_REAL *cout);

/* --------- other stuff -------- */
LTFAT_EXTERN void
LTFAT_NAME(heapint)(const LTFAT_REAL *s,
                    const LTFAT_REAL *tgrad,
                    const LTFAT_REAL *fgrad,
                    const ltfatInt a, const ltfatInt M,
                    const ltfatInt L, const ltfatInt W,
                    const LTFAT_REAL tol, LTFAT_REAL *phase);

LTFAT_EXTERN void
LTFAT_NAME(gabreassign)(const LTFAT_REAL *s, const LTFAT_REAL *tgrad,
                        const LTFAT_REAL *fgrad, const ltfatInt L, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M, LTFAT_REAL *sr);

LTFAT_EXTERN void
LTFAT_NAME(filterbankreassign)(const LTFAT_REAL*     s[],
                               const LTFAT_REAL* tgrad[],
                               const LTFAT_REAL* fgrad[],
                               const ltfatInt        N[],
                               const double          a[],
                               const double      cfreq[],
                               const ltfatInt          M,
                               LTFAT_REAL*          sr[],
                               fbreassHints        hints,
                               fbreassOptOut*      repos);

LTFAT_EXTERN void
LTFAT_NAME(filterbankphasegrad)(const LTFAT_COMPLEX* c [],
                                const LTFAT_COMPLEX* ch[],
                                const LTFAT_COMPLEX* cd[],
                                const ltfatInt          M,
                                const ltfatInt        N[],
                                const ltfatInt          L,
                                const LTFAT_REAL   minlvl,
                                LTFAT_REAL*        tgrad[],
                                LTFAT_REAL*        fgrad[],
                                LTFAT_REAL*           cs[]);

LTFAT_EXTERN void
LTFAT_NAME(fftshift_r)(const LTFAT_REAL *f, const ltfatInt L, LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(ifftshift_r)(const LTFAT_REAL *f, const ltfatInt L, LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(fir2long_r)(const LTFAT_REAL *f, const ltfatInt Lfir, const ltfatInt Llong,
                       LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(fir2long_c)(const LTFAT_COMPLEX *f,
                       const ltfatInt Lfir, const ltfatInt Llong,
                       LTFAT_COMPLEX *h);

LTFAT_EXTERN void
LTFAT_NAME(long2fir_r)(const LTFAT_REAL *f, const ltfatInt Llong,
                       const ltfatInt Lfir, LTFAT_REAL *h);

LTFAT_EXTERN void
LTFAT_NAME(long2fir_c)(const LTFAT_COMPLEX *f, const ltfatInt Llong,
                       const ltfatInt Lfir, LTFAT_COMPLEX *h);

ltfatInt
LTFAT_NAME(complexprod)(LTFAT_COMPLEX *c, const LTFAT_COMPLEX a,
                        const LTFAT_COMPLEX b);

/* ----- internal routines for calling BLAS and LAPACK ----- */

/*
// LAPACK overwrites the input argument.
ltfatInt
LTFAT_NAME(ltfat_posv)(const ltfatInt N, const ltfatInt NRHS,
             LTFAT_COMPLEX *A, const ltfatInt lda,
             LTFAT_COMPLEX *B, const ltfatInt ldb);

// LAPACK overwrites the input argument.
ltfatInt
LTFAT_NAME(ltfat_gesvd)(const ltfatInt M, const ltfatInt N,
              LTFAT_COMPLEX *A, const ltfatInt lda,
              LTFAT_REAL *S, LTFAT_COMPLEX *U, const ltfatInt ldu,
              LTFAT_COMPLEX *VT, const ltfatInt ldvt);

void
LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
             const enum CBLAS_TRANSPOSE TransB,
             const ltfatInt M, const ltfatInt N, const ltfatInt K,
             const LTFAT_COMPLEX *alpha,
             const LTFAT_COMPLEX *A, const ltfatInt lda,
             const LTFAT_COMPLEX *B, const ltfatInt ldb,
             const LTFAT_COMPLEX *beta,
             LTFAT_COMPLEX *C, const ltfatInt ldc);
*/


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

/*   --- dgtreal_long class definition  --- */
typedef struct
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt W;
    ltfatInt c;
    ltfatInt h_a;
    dgt_phasetype ptype;
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
                              const ltfatInt L, const ltfatInt W, const ltfatInt a,
                              const ltfatInt M, LTFAT_COMPLEX *cout,
                              const dgt_phasetype ptype, unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long_execute)(const LTFAT_NAME(dgtreal_long_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan) plan);




/*   --- dgt_fb class definition  --- */
typedef struct
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt gl;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_small;
    LTFAT_REAL *sbuf;
    LTFAT_REAL *fw;
    LTFAT_COMPLEX *gw;
    LTFAT_COMPLEX *cout;
} LTFAT_NAME(dgt_fb_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_fb_plan)
LTFAT_NAME(dgt_fb_init)(const LTFAT_COMPLEX *g,
                        const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_execute)(const LTFAT_NAME(dgt_fb_plan) plan,
                           const LTFAT_COMPLEX *f, const ltfatInt L,
                           const ltfatInt W, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_fb_done)(LTFAT_NAME(dgt_fb_plan) plan);



/*   --- dgtreal_fb class definition  --- */

typedef struct
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt gl;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_small;
    LTFAT_REAL    *sbuf;
    LTFAT_COMPLEX *cbuf;
    LTFAT_REAL *fw;
    LTFAT_REAL *gw;
    LTFAT_COMPLEX *cout;
} LTFAT_NAME(dgtreal_fb_plan);


LTFAT_EXTERN LTFAT_NAME(dgtreal_fb_plan)
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL *g,
                            const ltfatInt gl, const ltfatInt a,
                            const ltfatInt M, const dgt_phasetype ptype,
                            unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_execute)(const LTFAT_NAME(dgtreal_fb_plan) plan,
                               const LTFAT_REAL *f, const ltfatInt L, 
                               const ltfatInt W, LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan) plan);


/*   --- dgt_ola class definition  --- */
typedef struct
{
    LTFAT_NAME(dgt_long_plan) plan;
    ltfatInt bl;
    ltfatInt gl;
    ltfatInt W;
    LTFAT_COMPLEX *buf;
    LTFAT_COMPLEX *gext;
    LTFAT_COMPLEX *cbuf;

} LTFAT_NAME(dgt_ola_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_ola_plan)
LTFAT_NAME(dgt_ola_init)(const LTFAT_COMPLEX *g, const ltfatInt gl,
                         const ltfatInt W, const ltfatInt a,
                         const ltfatInt M, const ltfatInt bl,
                         const dgt_phasetype ptype, unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_execute)(const LTFAT_NAME(dgt_ola_plan) plan,
                            const LTFAT_COMPLEX *f, const ltfatInt L,
                            LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_done)(LTFAT_NAME(dgt_ola_plan) plan);


LTFAT_EXTERN void
LTFAT_NAME(dgt_walnut_plan)(LTFAT_NAME(dgt_long_plan) plan);


/*   --- dgtreal_ola class definition  --- */
typedef struct
{
    LTFAT_NAME(dgtreal_long_plan) plan;
    ltfatInt bl;
    ltfatInt gl;
    ltfatInt W;
    LTFAT_REAL *buf;
    LTFAT_REAL *gext;
    LTFAT_COMPLEX *cbuf;

} LTFAT_NAME(dgtreal_ola_plan);


LTFAT_EXTERN LTFAT_NAME(dgtreal_ola_plan)
LTFAT_NAME(dgtreal_ola_init)(const LTFAT_REAL *g, const ltfatInt gl,
                             const ltfatInt W, const ltfatInt a,
                             const ltfatInt M, const ltfatInt bl,
                             const dgt_phasetype ptype,
                             unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_execute)(const LTFAT_NAME(dgtreal_ola_plan) plan,
                                const LTFAT_REAL *f, const ltfatInt L,
                                LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_done)(LTFAT_NAME(dgtreal_ola_plan) plan);


LTFAT_EXTERN void
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan) plan);

/* -----  dgt_shearola class definition ------ */

typedef struct
{
    LTFAT_NAME(dgt_shear_plan) plan;
    ltfatInt bl;
    ltfatInt gl;
    ltfatInt W;
    LTFAT_COMPLEX *buf;
    LTFAT_COMPLEX *gext;
    LTFAT_COMPLEX *cbuf;

} LTFAT_NAME(dgt_shearola_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_shearola_plan)
LTFAT_NAME(dgt_shearola_init)(const LTFAT_COMPLEX *g, const ltfatInt gl,
                              const ltfatInt W, const ltfatInt a, const ltfatInt M,
                              const ltfatInt s0, const ltfatInt s1, const ltfatInt br,
                              const ltfatInt bl,
                              unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_execute)(const LTFAT_NAME(dgt_shearola_plan) plan,
                                 const LTFAT_COMPLEX *f, const ltfatInt L,
                                 LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_done)(LTFAT_NAME(dgt_shearola_plan) plan);



