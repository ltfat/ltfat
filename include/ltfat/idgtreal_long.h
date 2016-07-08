typedef struct LTFAT_NAME(idgtreal_long_plan) LTFAT_NAME(idgtreal_long_plan);

/** \defgroup idgtreallong Inverse Discrete Gabor Transform for real signals and long windows -- idgtreal_long
 *  \addtogroup idgtreallong
 * @{
 *
 *
 */


LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M,
                          const dgt_phasetype ptype, LTFAT_REAL *f);

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_init)(const LTFAT_COMPLEX *cin, const LTFAT_REAL *g,
                               const ltfatInt L, const ltfatInt W,
                               const ltfatInt a, const ltfatInt M, LTFAT_REAL *f,
                               const dgt_phasetype ptype, unsigned flags,
                               LTFAT_NAME(idgtreal_long_plan)** pout);

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_execute)(LTFAT_NAME(idgtreal_long_plan)* plan);

LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_done)(LTFAT_NAME(idgtreal_long_plan)** plan);

/** @}*/

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_walnut_execute)(LTFAT_NAME(idgtreal_long_plan)* p);
