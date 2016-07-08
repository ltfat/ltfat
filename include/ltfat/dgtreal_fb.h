typedef struct LTFAT_NAME(dgtreal_fb_plan) LTFAT_NAME(dgtreal_fb_plan);

/** \defgroup dgtrealfb Discrete Gabor Transform for real signals and FIR windows -- dgtreal_fb 
 *  \addtogroup dgtrealfb
 * @{
 *
 * \anchor dgt
 *  \f[
 *  c(m,n) 
 *   = \sum_{l=0}^{L-1}\! f(l)
 *   \overline{g(l-na)} \me^{-\mi 2\pi l m/M } \,
 *  \f]
 *
 *
 *  \f[
 *  c(m,n) 
 *   = \sum_{l=0}^{L-1}\! f(l)
 *   \overline{g(l-na)} \me^{-\mi 2\pi (l-na) m/M } \,
 *  \f]
 *
 */


LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                       const ltfatInt L, const ltfatInt gl,
                       const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                       const dgt_phasetype ptype, LTFAT_COMPLEX *cout);


LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL *g,
                            const ltfatInt gl, const ltfatInt a,
                            const ltfatInt M, const dgt_phasetype ptype,
                            unsigned flags, LTFAT_NAME(dgtreal_fb_plan)** pout);

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_execute)(LTFAT_NAME(dgtreal_fb_plan)* plan,
                               const LTFAT_REAL *f, const ltfatInt L,
                               const ltfatInt W, LTFAT_COMPLEX *cout);

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan)** plan);

/** @}*/
