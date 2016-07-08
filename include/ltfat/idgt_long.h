typedef struct LTFAT_NAME(idgt_long_plan) LTFAT_NAME(idgt_long_plan);

/** \defgroup idgtlong Inverse Discrete Gabor Transform for long windows -- idgt_long 
 *  \addtogroup idgtlong
 * @{
 *
 *
 */


LTFAT_EXTERN int
LTFAT_NAME(idgt_long)(const LTFAT_COMPLEX *cin, const LTFAT_TYPE *g,
                      const ltfatInt L, const ltfatInt W,
                      const ltfatInt a, const ltfatInt M,
                      const dgt_phasetype ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgt_fac)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *gf,
                     const ltfatInt L,
                     const ltfatInt W, const ltfatInt a, const ltfatInt M,
                     const dgt_phasetype ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN int
LTFAT_NAME(idgt_long_init)(const LTFAT_COMPLEX* cin, const LTFAT_TYPE* g,
                           const ltfatInt L, const ltfatInt W,
                           const ltfatInt a, const ltfatInt M,  LTFAT_COMPLEX* f,
                           unsigned flags, const dgt_phasetype ptype, LTFAT_NAME(idgt_long_plan)** pout);

LTFAT_EXTERN int
LTFAT_NAME(idgt_long_execute)(LTFAT_NAME(idgt_long_plan)* p);

LTFAT_EXTERN int
LTFAT_NAME(idgt_long_done)(LTFAT_NAME(idgt_long_plan)** plan);


/** @}*/

LTFAT_EXTERN void
LTFAT_NAME(idgt_walnut_execute)(LTFAT_NAME(idgt_long_plan)* p);
