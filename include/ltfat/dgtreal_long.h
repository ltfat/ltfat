/*   --- dgtreal_long class definition  --- */
typedef struct LTFAT_NAME(dgtreal_long_plan) LTFAT_NAME(dgtreal_long_plan);

/** \defgroup dgtreallong Discrete Gabor Transform for real signals and long windows -- dgtreal_long 
 *  \addtogroup dgtreallong
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
LTFAT_NAME(dgtreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                         const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                         const ltfatInt M, const dgt_phasetype ptype,
                         LTFAT_COMPLEX *cout);

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_init)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                              const ltfatInt L, const ltfatInt W, const ltfatInt a,
                              const ltfatInt M, LTFAT_COMPLEX *cout,
                              const dgt_phasetype ptype, unsigned flags,
                              LTFAT_NAME(dgtreal_long_plan)** plan);

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_execute)(LTFAT_NAME(dgtreal_long_plan)* plan);

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan)** plan);

/** @}*/

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan)* plan);
