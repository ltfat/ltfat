typedef struct LTFAT_NAME(idgtreal_fb_plan) LTFAT_NAME(idgtreal_fb_plan);

/** \defgroup idgtrealfb Inverse Discrete Gabor Transform for real signals and FIR windows -- idgtreal_fb 
 *  \addtogroup idgtrealfb
 * @{
 *
 *
 */

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                        const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, LTFAT_REAL *f);

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_init)(const LTFAT_REAL *g, const ltfatInt gl,
                             const ltfatInt a, const ltfatInt M, const dgt_phasetype ptype,
                             unsigned flags, LTFAT_NAME(idgtreal_fb_plan)** pout);

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_execute)(LTFAT_NAME(idgtreal_fb_plan)* p, const LTFAT_COMPLEX *cin,
                                const ltfatInt L, const ltfatInt W, LTFAT_REAL *f);

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_done)(LTFAT_NAME(idgtreal_fb_plan)** p);

/** @}*/
