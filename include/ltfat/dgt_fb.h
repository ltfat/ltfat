typedef struct LTFAT_NAME(dgt_fb_plan) LTFAT_NAME(dgt_fb_plan);

/** \defgroup dgtfb Discrete Gabor Transform for FIR windows -- dgt_fb 
 *  \addtogroup dgtfb
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
LTFAT_NAME(dgt_fb)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
        const ltfatInt L, const ltfatInt gl,
        const ltfatInt W,  const ltfatInt a, const ltfatInt M,
        const dgt_phasetype ptype, LTFAT_COMPLEX *cout);

LTFAT_EXTERN int
LTFAT_NAME(dgt_fb_init)(const LTFAT_TYPE *g,
                        const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, unsigned flags, LTFAT_NAME(dgt_fb_plan)** p);

LTFAT_EXTERN int
LTFAT_NAME(dgt_fb_execute)(const LTFAT_NAME(dgt_fb_plan)* plan,
                           const LTFAT_TYPE *f, const ltfatInt L,
                           const ltfatInt W, LTFAT_COMPLEX *cout);

LTFAT_EXTERN int
LTFAT_NAME(dgt_fb_done)(LTFAT_NAME(dgt_fb_plan)** plan);



/** @}*/
