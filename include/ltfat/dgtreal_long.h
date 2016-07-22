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

/** Compute Discrete Gabor Transform for real signals using the factorization algorithm
 *
 * \param[in]     f   Input signal, size L x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[out]    c   DGT coefficients, size M2 x N x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_long_d(const double f[], const double g[],
 *                      const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                      const ltfatInt M, const ltfat_phaseconvention ptype, complex double c[]);
 *
 * ltfat_dgtreal_long_s(const float f[], const float g[],
 *                      const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                      const ltfatInt M, const ltfat_phaseconvention ptype, complex float c[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long)(const LTFAT_REAL f[], const LTFAT_REAL g[],
                         const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                         const ltfatInt M, const ltfat_phaseconvention ptype,
                         LTFAT_COMPLEX c[]);

/** Initialize plan for Discrete Gabor Transform for real signals for the factorization algorithm
 *
 * \param[in]     f   Input signal, size L x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in]     c   DGT coefficients, size M2 x N x W
 * \param[in] ptype   Phase convention
 * \param[in] flags   FFTW plan flags
 * \param[out] plan   DGT plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_long_init_d(const double f[], const double g[],
 *                           const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                           const ltfatInt M,  complex double c[], const ltfat_phaseconvention ptype,
 *                           unsigned flags, ltfat_dgtreal_long_plan_d** plan);
 *
 * ltfat_dgtreal_long_init_s(const float f[], const float g[],
 *                           const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                           const ltfatInt M,  complex double c[], const ltfat_phaseconvention ptype,
 *                           unsigned flags, ltfat_dgtreal_long_plan_s** plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_init)(const LTFAT_REAL f[], const LTFAT_REAL g[],
                              const ltfatInt L, const ltfatInt W, const ltfatInt a,
                              const ltfatInt M, LTFAT_COMPLEX c[],
                              const ltfat_phaseconvention ptype, unsigned flags,
                              LTFAT_NAME(dgtreal_long_plan)** plan);

/** Execute plan for Discrete Gabor Transform for real signals using the factorization algorithm
 *
 * \param[in]  plan   DGT plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_long_execute_d(ltfat_dgtreal_long_plan_d* plan);
 *
 * ltfat_dgtreal_long_execute_s(ltfat_dgtreal_long_plan_s* plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_execute)(LTFAT_NAME(dgtreal_long_plan)* plan);

/** Execute plan for Discrete Gabor Transform for real signals using the factorization algorithm
 * 
 * ... on arrays which might have not been used in the init function.
 *
 * \param[in]  plan   DGT plan
 * \param[in]     f   Input signal, size L x W
 * \param[out]    c   Coefficients, size M2 x N x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_long_execute_newarray_d(ltfat_dgtreal_long_plan_d* plan, 
 *                                       const double f[], complex double c[]);
 *
 * ltfat_dgtreal_long_execute_newarray_s(ltfat_dgtreal_long_plan_s* plan
 *                                       const float f[], complex float c[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_execute_newarray)(LTFAT_NAME(dgtreal_long_plan)* plan,
                                          const LTFAT_REAL* f, LTFAT_COMPLEX* c);

/** Destroy the plan
 *
 * \param[in]  plan   DGT plan
 * \return Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_long_done_d(ltfat_dgtreal_long_plan_d** plan);
 *
 * ltfat_dgtreal_long_done_s(ltfat_dgtreal_long_plan_s** plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_long_done)(LTFAT_NAME(dgtreal_long_plan)** plan);

/** @}*/

LTFAT_EXTERN int
LTFAT_NAME(dgtreal_walnut_plan)(LTFAT_NAME(dgtreal_long_plan)* plan);
