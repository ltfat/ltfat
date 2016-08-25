#include "ltfat/types.h"

#include "dgt_long.h"
#include "idgt_long.h"
#include "dgt_fb.h"
#include "idgt_fb.h"
#include "wavelets.h"
#include "goertzel.h"
#include "ciutils.h"
#include "gabdual_painless.h"
#include "ci_windows.h"
#include "dct.h"
#include "dst.h"
#include "ci_memalloc.h"

/*   Walnut factorization    */

typedef struct LTFAT_NAME(wfac_plan) LTFAT_NAME(wfac_plan);

LTFAT_API int
LTFAT_NAME(wfac)(const LTFAT_TYPE *g, ltfat_int L, ltfat_int R,
                 ltfat_int a, ltfat_int M, LTFAT_COMPLEX *gf);

LTFAT_API int
LTFAT_NAME(wfac_init)(ltfat_int L, ltfat_int a, ltfat_int M,
                      unsigned flags, LTFAT_NAME(wfac_plan)** plan);

LTFAT_API int
LTFAT_NAME(wfac_execute)(LTFAT_NAME(wfac_plan)* plan, const LTFAT_TYPE *g,
                         ltfat_int R, LTFAT_COMPLEX *gf);

LTFAT_API int
LTFAT_NAME(wfac_done)(LTFAT_NAME(wfac_plan)** plan);

/*  Inverse Walnut factorization  */

typedef struct LTFAT_NAME(iwfac_plan) LTFAT_NAME(iwfac_plan);

LTFAT_API int
LTFAT_NAME(iwfac)(const LTFAT_COMPLEX *gf, ltfat_int L, ltfat_int R,
                  ltfat_int a, ltfat_int M, LTFAT_TYPE *g);

LTFAT_API int
LTFAT_NAME(iwfac_init)(ltfat_int L, ltfat_int a, ltfat_int M,
                       unsigned flags, LTFAT_NAME(iwfac_plan)** plan);

LTFAT_API int
LTFAT_NAME(iwfac_execute)(LTFAT_NAME(iwfac_plan)* plan, const LTFAT_COMPLEX* gf,
                          ltfat_int R, LTFAT_TYPE* g);

LTFAT_API int
LTFAT_NAME(iwfac_done)(LTFAT_NAME(iwfac_plan)** plan);


LTFAT_API void
LTFAT_NAME(col2diag)(const LTFAT_TYPE *cin, ltfat_int L,
                     LTFAT_TYPE *cout);

/*  Dual and tight  */

/** \addtogroup gabdual
 * @{
 *
 * In order to be able to compute canonical dual or tight windows,
 * the system must be a frame i.e. M>=a && gl>=a where gl is length
 * of the window support.
 *
 * In general, the canonical dual (tight) window support might not be the same
 * as of the original window.
 *
 * The canonical dual (tight) window is guaranteed to have the same support
 * in the _painless case_ i.e. if M>=gl holds in addition to the frame condition.
 *
 */


/** Compute canonical dual window for Gabor system
 *
 * \warning This function is not available if libltfat has been compiled with
 * NOBLASLAPACK.
 *
 * \param[in]   g    Original window(s), size L x R
 * \param[in]   L    Length of the system
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gd    Canonical dual window(s), size L x R
 *
 * #### Versions #
 * <tt>
 * ltfat_gabdual_long_d(const double g[], const ltfat_int L, 
 *                      ltfat_int a, ltfat_int M, double gd[]);
 *
 * ltfat_gabdual_long_s(const float g[], ltfat_int L, 
 *                      ltfat_int a, ltfat_int M, float gd[]);
 *
 * ltfat_gabdual_long_dc(const ltfat_complex_d g[], ltfat_int L, 
 *                       ltfat_int a, ltfat_int M, ltfat_complex_d gd[]);
 *
 * ltfat_gabdual_long_sc(const ltfat_complex_s g[], ltfat_int L, 
 *                       ltfat_int a, ltfat_int M, ltfat_complex_s gd[]);
 * </tt>
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Either of the arrays is NULL
 * LTFATERR_BADSIZE         | Size of the array \a L is less or equal to 0.
 * LTFATERR_BADTRALEN       | \a L is not divisible by both \a a and \a M.
 * LTFATERR_NOTPOSARG       | Either of \a R, \a a, \a M is less or equal to 0.
 * LTFATERR_NOTAFRAME       | System does not form a frame i.e. M<a
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(gabdual_long)(const LTFAT_TYPE g[],
                         ltfat_int L, ltfat_int a,
                         ltfat_int M, LTFAT_TYPE gd[]);

/** Compute canonical tight window for Gabor system
 *
 * \warning This function is not available if libltfat has been compiled with
 * NOBLASLAPACK.
 * 
 * \see fir2long long2fir
 *
 * \param[in]   g    Original window(s), size L x R
 * \param[in]   L    Length of the system
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gt    Canonical tight window
 *
 * #### Versions #
 * ltfat_gabtight_long_d(const double g[], ltfat_int L, 
 *                       ltfat_int a, ltfat_int M, double gt[]);
 *
 * ltfat_gabtight_long_s(const float g[], ltfat_int L, 
 *                       ltfat_int a, ltfat_int M, float gt[]);
 *
 * ltfat_gabtight_long_dc(const ltfat_complex_d g[], ltfat_int L, 
 *                        ltfat_int a, ltfat_int M, ltfat_complex_d gt[]);
 *
 * ltfat_gabtight_long_sc(const ltfat_complex_s g[], ltfat_int L, 
 *                        ltfat_int a, ltfat_int M, ltfat_complex_s gt[]);
 * </tt>
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Either of the arrays is NULL
 * LTFATERR_BADSIZE         | Size of the array \a L is less or equal to 0.
 * LTFATERR_BADTRALEN       | L is not divisible by both \a a and \a M.
 * LTFATERR_NOTPOSARG       | Either of \a R, \a a, \a M is less or equal to 0.
 * LTFATERR_NOTAFRAME       | System does not form a frame i.e. M<a
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(gabtight_long)(const LTFAT_TYPE g[],
                          ltfat_int L, ltfat_int a,
                          ltfat_int M, LTFAT_TYPE gd[]);
/** @} */

/** Compute FIR canonical dual window for Gabor system
 *
 * The function internally calls fir2long, gabdual_long and long2fir. The window
 * might no longer be exact canonical dual window if gdl is smaller than the
 * length of the support of the window.
 *
 * \warning This function is not available if libltfat has been compiled with
 * NOBLASLAPACK.
 *
 * \param[in]    g    Original window
 * \param[in]   gl    Length of the window
 * \param[in]    L    Length of the system
 * \param[in]    a    Hop factor
 * \param[in]    M    Number of channels
 * \param[in]  gdl    Length of the dual window
 * \param[out]  gd    Canonical dual window
 *
 * #### Versions #
 * <tt>
 * ltfat_gabdual_fir_d(const double g[], ltfat_int gl, ltfat_int L,
 *                     ltfat_int a, ltfat_int M, ltfat_int gdl,
 *                     double gd[]);
 *
 * ltfat_gabdual_fir_s(const float g[], ltfat_int gl, ltfat_int L,
 *                     ltfat_int a, ltfat_int M, ltfat_int gdl,
 *                     float gd[]);
 *
 * ltfat_gabdual_fir_dc(const ltfat_complex_d g[], ltfat_int gl, ltfat_int L,
 *                      ltfat_int a, ltfat_int M, ltfat_int gdl,
 *                      ltfat_complex_d gd[]);
 *
 * ltfat_gabdual_fir_sc(const ltfat_complex_s g[], ltfat_int gl, ltfat_int L,
 *                      ltfat_int a, ltfat_int M, ltfat_int gdl,
 *                      ltfat_complex_s gd[]);
 * </tt>
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Either of the arrays is NULL
 * LTFATERR_BADSIZE         | Either of the array sizes is less or equal to 0.
 * LTFATERR_NOTPOSARG       | Either of \a a, \a M is less or equal to 0.
 * LTFATERR_BADTRALEN       | L is not divisible by both \a a and \a M.
 * LTFATERR_BADREQSIZE      | \a L is not greater or equal than both \a gl and \a gdl
 * LTFATERR_NOTAFRAME       | System does not form a frame i.e. M<a
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(gabdual_fir)(const LTFAT_TYPE g[], ltfat_int gl,
                        ltfat_int L, ltfat_int a,
                        ltfat_int M, ltfat_int gdl, LTFAT_TYPE gd[]);

/** Compute FIR canonical tight window for Gabor system
 *
 * The function internally calls fir2long, gabtight_long and long2fir. The window
 * might no longer be exact canonical tight window if gdl is smaller than the
 * length of the support of the window.
 *
 * \warning This function is not available if libltfat has been compiled with
 * NOBLASLAPACK.
 *
 * \param[in]    g    Original window, size gl x 1
 * \param[in]   gl    Length of the window
 * \param[in]    L    Length of the system
 * \param[in]    a    Hop factor
 * \param[in]    M    Number of channels
 * \param[in]  gtl    Length of the tight window
 * \param[out]  gt    Canonical dual window, size gtl x 1
 *
 * #### Versions #
 * <tt>
 * ltfat_gabtight_fir_d(const double g[], ltfat_int gl, ltfat_int L,
 *                      ltfat_int a, ltfat_int M, ltfat_int gtl,
 *                      double gt[]);
 *
 * ltfat_gabtight_fir_s(const float g[], ltfat_int gl, ltfat_int L,
 *                      ltfat_int a, ltfat_int M, ltfat_int gtl,
 *                      float gt[]);
 *
 * ltfat_gabtight_fir_dc(const ltfat_complex_d g[], ltfat_int gl, ltfat_int L,
 *                       ltfat_int a, ltfat_int M, ltfat_int gtl,
 *                       ltfat_complex_d gt[]);
 *
 * ltfat_gabtight_fir_sc(const ltfat_complex_s g[], ltfat_int gl, ltfat_int L,
 *                       ltfat_int a, ltfat_int M, ltfat_int gtl,
 *                       ltfat_complex_s gt[]);
 * </tt>
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Either of the arrays is NULL
 * LTFATERR_BADSIZE         | Either of the array sizes is less or equal to 0.
 * LTFATERR_NOTPOSARG       | Either of \a a, \a M is less or equal to 0.
 * LTFATERR_BADTRALEN       | L is not divisible by both \a a and \a M.
 * LTFATERR_BADREQSIZE      | \a L is not greater or equal than both \a gl and \a gtl
 * LTFATERR_NOTAFRAME       | System does not form a frame i.e. M<a
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(gabtight_fir)(const LTFAT_TYPE g[], ltfat_int gl,
                         ltfat_int L, ltfat_int a,
                         ltfat_int M, ltfat_int gtl, LTFAT_TYPE gt[]);


/* --------- Wilson and WMDCT bases ---------*/
LTFAT_API void
LTFAT_NAME(dwilt_long)(const LTFAT_TYPE *f,
                       const LTFAT_TYPE *g,
                       ltfat_int L, ltfat_int W, ltfat_int M,
                       LTFAT_TYPE *cout);

LTFAT_API void
LTFAT_NAME(dwilt_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                     ltfat_int L, ltfat_int gl, ltfat_int W, ltfat_int M,
                     LTFAT_TYPE *cout);


LTFAT_API void
LTFAT_NAME(dwiltiii_long)(const LTFAT_TYPE *f,
                          const LTFAT_TYPE *g,
                          ltfat_int L, ltfat_int W, ltfat_int M,
                          LTFAT_TYPE *cout);

LTFAT_API void
LTFAT_NAME(dwiltiii_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                        ltfat_int L, ltfat_int gl, ltfat_int W, ltfat_int M,
                        LTFAT_TYPE *cout);


/* --------- Wilson and WMDCT inverses ---------*/


LTFAT_API void
LTFAT_NAME(idwilt_long)(const LTFAT_TYPE *cin,
                        const LTFAT_TYPE *g,
                        ltfat_int L, ltfat_int W, ltfat_int M,
                        LTFAT_TYPE *f);

LTFAT_API void
LTFAT_NAME(idwilt_fb)(const LTFAT_TYPE *cin, const LTFAT_TYPE *g,
                      ltfat_int L, ltfat_int gl, ltfat_int W, ltfat_int M,
                      LTFAT_TYPE *f);

LTFAT_API void
LTFAT_NAME(idwiltiii_long)(const LTFAT_TYPE *cin,
                           const LTFAT_TYPE *g,
                           ltfat_int L, ltfat_int W, ltfat_int M,
                           LTFAT_TYPE *f);

LTFAT_API void
LTFAT_NAME(idwiltiii_fb)(const LTFAT_TYPE *cin, const LTFAT_TYPE *g,
                         ltfat_int L, ltfat_int gl, ltfat_int W, ltfat_int M,
                         LTFAT_TYPE *f);

/* --------------- DCT -------------------*/

// LTFAT_API LTFAT_FFTW(plan)
// LTFAT_NAME(dct_init)( ltfat_int L, ltfat_int W, LTFAT_TYPE *cout,
//                       const dct_kind kind);
//
//
// LTFAT_API void
// LTFAT_NAME(dct)(const LTFAT_TYPE *f, ltfat_int L, ltfat_int W,
//                 LTFAT_TYPE *cout, const dct_kind kind);
//
// LTFAT_API void
// LTFAT_NAME(dct_execute)(const LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
//                         ltfat_int L, ltfat_int W,
//                         LTFAT_TYPE *cout, const dct_kind kind);

/* --------------- DST -------------------*/

// LTFAT_API LTFAT_FFTW(plan)
// LTFAT_NAME(dst_init)( ltfat_int L, ltfat_int W, LTFAT_TYPE *cout,
//                       const dst_kind kind);
//
// LTFAT_API void
// LTFAT_NAME(dst)(const LTFAT_TYPE *f, ltfat_int L, ltfat_int W,
//                 LTFAT_TYPE *cout, const dst_kind kind);
//
// LTFAT_API void
// LTFAT_NAME(dst_execute)(LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
//                         ltfat_int L, ltfat_int W, LTFAT_TYPE *cout,
//                         const dst_kind kind);

/* --------------- Reassignment -----------*/

LTFAT_API void
LTFAT_NAME(gabreassign)(const LTFAT_TYPE *s, const LTFAT_REAL *tgrad,
                        const LTFAT_REAL *fgrad, ltfat_int L, ltfat_int W,
                        ltfat_int a, ltfat_int M, LTFAT_TYPE *sr);

LTFAT_API void
LTFAT_NAME(filterbankreassign)(const LTFAT_TYPE*     s[],
                               const LTFAT_REAL* tgrad[],
                               const LTFAT_REAL* fgrad[],
                               ltfat_int        N[],
                               const double          a[],
                               const double      cfreq[],
                               ltfat_int          M,
                               LTFAT_TYPE*          sr[],
                               fbreassHints        hints,
                               fbreassOptOut*      repos);
