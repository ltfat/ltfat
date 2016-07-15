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

/** Compute Discrete Gabor Transform for real signals using filter bank algorithm
 *
 * \param[in]     f   Input signal, size L x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]    gl   Window length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[out]    c   DGT coefficients, size M2 x N x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_fb_d(const double f[], const double g[],
 *                    const ltfatInt L, const ltfatInt gl,
 *                    const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                    const dgt_phasetype ptype, complex double c[]);
 *
 * ltfat_dgtreal_fb_s(const float f[], const float g[],
 *                    const ltfatInt L, const ltfatInt gl,
 *                    const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                    const dgt_phasetype ptype, complex float c[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb)(const LTFAT_REAL f[], const LTFAT_REAL g[],
                       const ltfatInt L, const ltfatInt gl,
                       const ltfatInt W,  const ltfatInt a, const ltfatInt M,
                       const dgt_phasetype ptype, LTFAT_COMPLEX c[]);

/** Initialize plan for Discrete Gabor Transform for real signals for the filter bank algorithm
 *
 * \param[in]     g   Window, size L x 1
 * \param[in]    gl   Window length
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[in] flags   FFTW plan flags
 * \param[out] pout   DGT plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_fb_init_d(const double g[], const ltfatInt gl, const ltfatInt a,
 *                         const ltfatInt M, const dgt_phasetype ptype, unsigned flags
 *                         ltfat_dgtreal_fb_plan_d** pout);
 *
 * ltfat_dgtreal_fb_init_s(const float g[], const ltfatInt gl, const ltfatInt a,
 *                         const ltfatInt M, const dgt_phasetype ptype,unsigned flags
 *                         ltfat_dgtreal_fb_plan_s** pout);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_init)(const LTFAT_REAL g[],
                            const ltfatInt gl, const ltfatInt a,
                            const ltfatInt M, const dgt_phasetype ptype,
                            unsigned flags, LTFAT_NAME(dgtreal_fb_plan)** pout);

/** Execute plan for Discrete Gabor Transform for real signals using the filter bank algorithm
 *
 * \param[in]  plan   DGT plan
 * \param[in]     f   Input signal, size L x W
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[out]    c   DGT coefficients, size M2 x N x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_fb_execute_d(ltfat_dgtreal_fb_plan_d* plan, const double f[],
 *                            const ltfatInt L, const ltfatInt W, complex double c[]);
 *
 * ltfat_dgtreal_fb_execute_s(ltfat_dgtreal_fb_plan_s* plan, const float f[],
 *                            const ltfatInt L, const ltfatInt W, complex float c[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_execute)(LTFAT_NAME(dgtreal_fb_plan)* plan,
                               const LTFAT_REAL f[], const ltfatInt L,
                               const ltfatInt W, LTFAT_COMPLEX c[]);

/** Destroy the plan
 *
 * \param[in]  plan   DGT plan
 * \return Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_dgtreal_fb_done_d(ltfat_dgtreal_fb_plan_d** plan);
 *
 * ltfat_dgtreal_fb_done_s(ltfat_dgtreal_fb_plan_s** plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgtreal_fb_done)(LTFAT_NAME(dgtreal_fb_plan)** plan);

/** @}*/
