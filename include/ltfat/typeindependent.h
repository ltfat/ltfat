#include "dgtreal_long.h"
#include "idgtreal_long.h"
#include "dgtreal_fb.h"
#include "idgtreal_fb.h"
#include "dgt_multi.h"
#include "dgt_shear.h"
#include "tiutils.h"
#include "rtdgtreal.h"
#include "fftw_wrappers.h"

/*  --------- factorizations --------------- */

LTFAT_EXTERN void
LTFAT_NAME(wfacreal)(const LTFAT_REAL *g, const ltfatInt L, const ltfatInt R,
                     const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *gf);

LTFAT_EXTERN void
LTFAT_NAME(iwfacreal)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                      const ltfatInt a, const ltfatInt M, LTFAT_REAL *g);

/* --------- DGT by factorization ------------ */

// LTFAT_EXTERN void
// LTFAT_NAME(dgt_fac)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *gf,
//                     const ltfatInt L, const ltfatInt W,  const ltfatInt a,
//                     const ltfatInt M, LTFAT_COMPLEX *cout);


// LTFAT_EXTERN void
// LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
//                          const ltfatInt L, const ltfatInt W,  const ltfatInt a,
//                          const ltfatInt M, const ltfat_phaseconvention ptype,
//                          LTFAT_COMPLEX *cout);

// LTFAT_EXTERN void
// LTFAT_NAME(dgt_fac_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
//                       const ltfatInt L, const ltfatInt W,  const ltfatInt a,
//                       const ltfatInt M, const ltfat_phaseconvention ptype,
//                       LTFAT_COMPLEX *cout);

// LTFAT_EXTERN
// void LTFAT_NAME(dgtreal_fac)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
//                              const ltfatInt L,
//                              const ltfatInt W,  const ltfatInt a,
//                              const ltfatInt M, LTFAT_COMPLEX *cout);

// LTFAT_EXTERN void
// LTFAT_NAME(dgt_walnut_r)(const LTFAT_REAL *f, const LTFAT_COMPLEX *gf,
//                          const ltfatInt L, const ltfatInt W,
//                          const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *cout);





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


// LTFAT_EXTERN void
// LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
//                        const ltfatInt L, const ltfatInt gl,
//                        const ltfatInt W,  const ltfatInt a, const ltfatInt M,
//                        const ltfat_phaseconvention ptype, LTFAT_COMPLEX *cout);

// LTFAT_EXTERN void
// LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_COMPLEX *g,
//                     const ltfatInt L, const ltfatInt gl,
//                     const ltfatInt W, const ltfatInt a, const ltfatInt M,
//                     const ltfat_phaseconvention ptype, LTFAT_COMPLEX *f);
//
// LTFAT_EXTERN void
// LTFAT_NAME(idgt_fb_r)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
//                       const ltfatInt L, const ltfatInt gl,
//                       const ltfatInt W, const ltfatInt a, const ltfatInt M,
//                       LTFAT_COMPLEX *f);

/* ---------- OLA DGTs ------------- */
LTFAT_EXTERN void
LTFAT_NAME(dgt_ola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a,
                    const ltfatInt M, const ltfatInt bl,
                    const ltfat_phaseconvention ptype,
                    LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                        const ltfatInt L, const ltfatInt gl,
                        const ltfatInt W, const ltfatInt a, const ltfatInt M,
                        const ltfatInt bl, const ltfat_phaseconvention ptype,
                        LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                         const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                         const ltfatInt a, const ltfatInt M,
                         const ltfatInt s0, const ltfatInt s1, const ltfatInt br,
                         const ltfatInt bl, LTFAT_COMPLEX *cout);

/* --------- FFT ------------------*/
// LTFAT_EXTERN LTFAT_FFTW(plan)
// LTFAT_NAME(fftreal_init)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
//                          LTFAT_COMPLEX *cout, unsigned flag);
//
// LTFAT_EXTERN void
// LTFAT_NAME(fftreal_execute)(LTFAT_FFTW(plan) p, LTFAT_REAL *f,
//                             LTFAT_COMPLEX *cout);
//
// LTFAT_EXTERN void
// LTFAT_NAME(fftreal)(LTFAT_REAL *f, const ltfatInt L, const ltfatInt W,
//                     LTFAT_COMPLEX *cout);
//
// LTFAT_EXTERN LTFAT_FFTW(plan)
// LTFAT_NAME(ifftreal_init)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
//                           LTFAT_REAL *f, unsigned flag);
//
// LTFAT_EXTERN void
// LTFAT_NAME(ifftreal_execute)(LTFAT_FFTW(plan), LTFAT_COMPLEX *c,
//                              const ltfatInt L, const ltfatInt W,
//                              LTFAT_REAL *f);
//
// LTFAT_EXTERN void
// LTFAT_NAME(ifftreal)(LTFAT_COMPLEX *c, const ltfatInt L, const ltfatInt W,
//                      LTFAT_REAL *f);


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
/** \addtogroup windows
 * @{
 */

/** Compute real, periodized Gaussian window
 *
 * \param[in]   L      Window length
 * \param[in]   w      Time-freqency support ratio
 * \param[in]   c_t    Time center offset
 * \param[out]  g      Window
 *
 * #### Function versions #
 * <tt>
 * ltfat_pgauss_d(const ltfatInt L,const double w, const double c_t, double* g);
 *
 * ltfat_pgauss_s(const ltfatInt L,const double w, const double c_t, float* g);
 * </tt>
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | The output array is NULL.
 * LTFATERR_BADSIZE         | Window length is less or equal to 0.
 * LTFATERR_NOTPOSARG    | \a w is less or equal to zero.
 */
LTFAT_EXTERN int
LTFAT_NAME(pgauss)(const ltfatInt L, const double w, const double c_t,
                   LTFAT_REAL *g);

/** Compute complex, periodized Gaussian window
 *
 * \param[in]   L      Window length
 * \param[in]   w      Time-freqency support ratio
 * \param[in]   c_t    Time center offset
 * \param[in]   c_f    Frequency center offset
 * \param[out]  g      Window
 *
 * #### Function versions #
 *
 * <tt>
 * ltfat_pgauss_cmplx_d(const ltfatInt L, const double w, const double c_t,
 *                      const double c_f, complex double* g);
 *
 * ltfat_pgauss_cmplx_s(const ltfatInt L, const double w, const double c_t,
 *                      const double c_f, complex float* g);
 * </tt>
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | The output array is NULL.
 * LTFATERR_BADSIZE         | Window length is less or equal to 0.
 * LTFATERR_NOTPOSARG    | \a w is less or equal to zero.
 */
LTFAT_EXTERN int
LTFAT_NAME(pgauss_cmplx)(const ltfatInt L, const double w, const double c_t,
                         const double c_f, LTFAT_COMPLEX *g);

/** @} */

/* --------- pfilt and filterbanks ------------- */
LTFAT_EXTERN void
LTFAT_NAME(pfilt_fir_rr)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                         const ltfatInt L, const ltfatInt gl,
                         const ltfatInt W, const ltfatInt a,
                         LTFAT_REAL *cout);


/* --------- other stuff -------- */

#ifndef _ltfat_mask_element_defined
#define _ltfat_mask_element_defined

enum ltfat_mask_element
{
    LTFAT_MASK_BELOWTOL    = -1, // Do not compute phase, the coefficient is too small
    LTFAT_MASK_UNKNOWN     =  0, // Will compute phase for these
    LTFAT_MASK_KNOWN       =  1, // The phase was already known
    LTFAT_MASK_WENTNORTH   =  2, // Phase was spread from the south neighbor
    LTFAT_MASK_WENTSOUTH   =  3, // Phase was spread from the north neighbor
    LTFAT_MASK_WENTEAST    =  4, // Phase was spread from the west neighbor
    LTFAT_MASK_WENTWEST    =  5, // Phase was spread from the east neighbor
    LTFAT_MASK_STARTPOINT  =  6, // This is the initial point of integration. It gets zero phase
    LTFAT_MASK_BORDERPOINT =  7, // This is candidate border coefficient with known phase
};

#endif

typedef struct LTFAT_NAME(heapinttask) LTFAT_NAME(heapinttask);

LTFAT_EXTERN LTFAT_NAME(heapinttask)*
LTFAT_NAME(heapinttask_init)(const ltfatInt height, const ltfatInt N,
                             const ltfatInt initheapsize,
                             const LTFAT_REAL* s, int do_real);

LTFAT_EXTERN void
LTFAT_NAME(heapint_execute)(LTFAT_NAME(heapinttask)* hit,
                            const LTFAT_REAL* tgradw,
                            const LTFAT_REAL* fgradw,
                            LTFAT_REAL* phase);

LTFAT_EXTERN void
LTFAT_NAME(heapinttask_done)(LTFAT_NAME(heapinttask)* hit);

LTFAT_EXTERN void
LTFAT_NAME(heapinttask_resetmax)(LTFAT_NAME(heapinttask)* hit,
                                 const LTFAT_REAL* news,
                                 const LTFAT_REAL tol);


LTFAT_EXTERN void
LTFAT_NAME(heapinttask_resetmask)(LTFAT_NAME(heapinttask)* hit,
                                  const int* mask,
                                  const LTFAT_REAL* news,
                                  const LTFAT_REAL tol,
                                  const int do_log);

LTFAT_EXTERN int*
LTFAT_NAME(heapinttask_get_mask)( LTFAT_NAME(heapinttask)* hit);

LTFAT_EXTERN void
LTFAT_NAME(heapint)(const LTFAT_REAL *s,
                    const LTFAT_REAL *tgradw,
                    const LTFAT_REAL *fgradw,
                    const ltfatInt a, const ltfatInt M,
                    const ltfatInt L, const ltfatInt W,
                    const LTFAT_REAL tol, LTFAT_REAL *phase);

// Does the same as the previous but
LTFAT_EXTERN void
LTFAT_NAME(heapint_relgrad)(const LTFAT_REAL *s,
                            const LTFAT_REAL *tgrad,
                            const LTFAT_REAL *fgrad,
                            const ltfatInt a, const ltfatInt M,
                            const ltfatInt L, const ltfatInt W,
                            const LTFAT_REAL tol, ltfat_phaseconvention phasetype,
                            LTFAT_REAL *phase);


LTFAT_EXTERN void
LTFAT_NAME(maskedheapint)(const LTFAT_REAL  *c,
                          const LTFAT_REAL *tgradw,
                          const LTFAT_REAL *fgradw,
                          const int* mask,
                          const ltfatInt a, const ltfatInt M,
                          const ltfatInt L, const ltfatInt W,
                          LTFAT_REAL tol, LTFAT_REAL *phase);

LTFAT_EXTERN void
LTFAT_NAME(maskedheapint_relgrad)(const LTFAT_REAL  *c,
                                  const LTFAT_REAL *tgrad,
                                  const LTFAT_REAL *fgrad,
                                  const int* mask,
                                  const ltfatInt a, const ltfatInt M,
                                  const ltfatInt L, const ltfatInt W,
                                  LTFAT_REAL tol, ltfat_phaseconvention phasetype,
                                  LTFAT_REAL *phase);

LTFAT_EXTERN void
LTFAT_NAME(heapintreal)(const LTFAT_REAL *s,
                        const LTFAT_REAL *tgradw,
                        const LTFAT_REAL *fgradw,
                        const ltfatInt a, const ltfatInt M,
                        const ltfatInt L, const ltfatInt W,
                        const LTFAT_REAL tol,
                        LTFAT_REAL *phase);

LTFAT_EXTERN void
LTFAT_NAME(heapintreal_relgrad)(const LTFAT_REAL *s,
                                const LTFAT_REAL *tgradw,
                                const LTFAT_REAL *fgradw,
                                const ltfatInt a, const ltfatInt M,
                                const ltfatInt L, const ltfatInt W,
                                const LTFAT_REAL tol, ltfat_phaseconvention phasetype,
                                LTFAT_REAL *phase);

LTFAT_EXTERN void
LTFAT_NAME(maskedheapintreal)(const LTFAT_REAL * s,
                              const LTFAT_REAL * tgrad,
                              const LTFAT_REAL * fgrad,
                              const int* mask,
                              const ltfatInt a, const ltfatInt M,
                              const ltfatInt L, const ltfatInt W,
                              LTFAT_REAL tol,
                              LTFAT_REAL * phase);

LTFAT_EXTERN void
LTFAT_NAME(maskedheapintreal_relgrad)(const LTFAT_REAL* s,
                                      const LTFAT_REAL* tgradw,
                                      const LTFAT_REAL* fgradw,
                                      const int* mask,
                                      const ltfatInt a, const ltfatInt M,
                                      const ltfatInt L, const ltfatInt W,
                                      LTFAT_REAL tol, ltfat_phaseconvention phasetype,
                                      LTFAT_REAL* phase);

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


// // LAPACK overwrites the input argument.
// ltfatInt
// LTFAT_NAME(ltfat_posv)(const ptrdiff_t N, const ptrdiff_t NRHS,
//                        LTFAT_COMPLEX *A, const ptrdiff_t lda,
//                        LTFAT_COMPLEX *B, const ptrdiff_t ldb);
//
// // LAPACK overwrites the input argument.
// ltfatInt
// LTFAT_NAME(ltfat_gesvd)(const ptrdiff_t M, const ptrdiff_t N,
//                         LTFAT_COMPLEX *A, const ptrdiff_t lda,
//                         LTFAT_REAL *S, LTFAT_COMPLEX *U, const ptrdiff_t ldu,
//                         LTFAT_COMPLEX *VT, const ptrdiff_t ldvt);
//
// void
// LTFAT_NAME(ltfat_gemm)(const enum CBLAS_TRANSPOSE TransA,
//                        const enum CBLAS_TRANSPOSE TransB,
//                        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
//                        const LTFAT_COMPLEX *alpha,
//                        const LTFAT_COMPLEX *A, const ptrdiff_t lda,
//                        const LTFAT_COMPLEX *B, const ptrdiff_t ldb,
//                        const LTFAT_COMPLEX *beta,
//                        LTFAT_COMPLEX *C, const ptrdiff_t ldc);


/*   --- dgtreal_fb class definition  --- */

// typedef struct
// {
//     ltfatInt a;
//     ltfatInt M;
//     ltfatInt gl;
//     ltfat_phaseconvention ptype;
//     LTFAT_FFTW(plan) p_small;
//     LTFAT_REAL    *sbuf;
//     LTFAT_COMPLEX *cbuf;
//     LTFAT_REAL *fw;
//     LTFAT_REAL *gw;
//     LTFAT_COMPLEX *cout;
// } LTFAT_NAME(dgtreal_fb_plan);
//
//
// LTFAT_EXTERN LTFAT_NAME(dgtreal_fb_plan)
// LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL *g,
//                             const ltfatInt gl, const ltfatInt a,
//                             const ltfatInt M, const ltfat_phaseconvention ptype,
//                             unsigned flags);
//
// LTFAT_EXTERN void
// LTFAT_NAME(dgtreal_fb_execute)(const LTFAT_NAME(dgtreal_fb_plan) plan,
//                                const LTFAT_REAL *f, const ltfatInt L,
//                                const ltfatInt W, LTFAT_COMPLEX *cout);
//
// LTFAT_EXTERN void
// LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan) plan);


/*   --- dgt_ola class definition  --- */
typedef struct
{
    LTFAT_NAME_COMPLEX(dgt_long_plan)* plan;
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
                         const ltfat_phaseconvention ptype, unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_execute)(const LTFAT_NAME(dgt_ola_plan) plan,
                            const LTFAT_COMPLEX *f, const ltfatInt L,
                            LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_ola_done)(LTFAT_NAME(dgt_ola_plan) plan);


// LTFAT_EXTERN void
// LTFAT_NAME(dgt_walnut_plan)(LTFAT_NAME(dgt_long_plan)* plan);


/*   --- dgtreal_ola class definition  --- */
typedef struct
{
    LTFAT_NAME(dgtreal_long_plan)* plan;
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
                             const ltfat_phaseconvention ptype,
                             unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_execute)(const LTFAT_NAME(dgtreal_ola_plan) plan,
                                const LTFAT_REAL *f, const ltfatInt L,
                                LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtreal_ola_done)(LTFAT_NAME(dgtreal_ola_plan) plan);


// LTFAT_EXTERN void
// LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan) plan);

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
