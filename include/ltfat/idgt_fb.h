typedef struct LTFAT_NAME(idgt_fb_plan) LTFAT_NAME(idgt_fb_plan);

/** 
 *  \addtogroup dgtfb
 * @{
 *
 *
 */

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX *cin, const LTFAT_TYPE *g,
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
                    const ltfat_phaseconvention ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb_init)(const LTFAT_TYPE *g, const ltfatInt gl,
                         const ltfatInt a, const ltfatInt M, const ltfat_phaseconvention ptype,
                         unsigned flags, LTFAT_NAME(idgt_fb_plan)** pout);

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb_execute)(LTFAT_NAME(idgt_fb_plan)* p, const LTFAT_COMPLEX *cin,
                            const ltfatInt L, const ltfatInt W, LTFAT_COMPLEX *f);

LTFAT_EXTERN int
LTFAT_NAME(idgt_fb_done)(LTFAT_NAME(idgt_fb_plan)** p);


/** @}*/
