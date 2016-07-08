typedef struct LTFAT_NAME(idgt_fb_plan) LTFAT_NAME(idgt_fb_plan);

/** \defgroup idgtfb Inverse Discrete Gabor Transform for FIR windows -- idgt_fb 
 *  \addtogroup idgtfb
 * @{
 *
 *
 */

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_TYPE *g,
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
                    const dgt_phasetype ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb_init)(const LTFAT_TYPE *g, const ltfatInt gl,
                         const ltfatInt a, const ltfatInt M, const dgt_phasetype ptype,
                         unsigned flags, LTFAT_NAME(idgt_fb_plan)** pout);

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb_execute)(LTFAT_NAME(idgt_fb_plan)* p, const LTFAT_COMPLEX *cin,
                            const ltfatInt L, const ltfatInt W, LTFAT_COMPLEX *f);

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb_done)(LTFAT_NAME(idgt_fb_plan)** p);


/** @}*/
