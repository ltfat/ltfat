typedef struct LTFAT_NAME(idgtreal_fb_plan) LTFAT_NAME(idgtreal_fb_plan);

/**
 *  \addtogroup dgtrealfb
 * @{
 *
 *
 */

/** Compute Inverse Discrete Gabor Transform for real signals using filter bank algorithm
 *
 * \param[out]    c   DGT coefficients, size M2 x N x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]    gl   Window length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[in]     f   Output signal, size L x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_fb_d(const complex double c[], const double g[],
 *                     const ltfatInt L, const ltfatInt gl,
 *                     const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                     const dgt_phasetype ptype, double f[]);
 *
 * ltfat_idgtreal_fb_s(const complex float c[], const float g[],
 *                     const ltfatInt L, const ltfatInt gl,
 *                     const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                     const dgt_phasetype ptype, float f[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb)(const LTFAT_COMPLEX c[], const LTFAT_REAL g[],
                        const ltfatInt L, const ltfatInt gl, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M,
                        const dgt_phasetype ptype, LTFAT_REAL f[]);

/** Initialize plan for Inverse Discrete Gabor Transform for real signals for the filter bank algorithm
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
 * ltfat_idgtreal_fb_init_d(const double g[], const ltfatInt gl, const ltfatInt a,
 *                          const ltfatInt M, const dgt_phasetype ptype, unsigned flags
 *                          ltfat_idgtreal_fb_plan_d** pout);
 *
 * ltfat_idgtreal_fb_init_s(const float g[], const ltfatInt gl, const ltfatInt a,
 *                          const ltfatInt M, const dgt_phasetype ptype, unsigned flags
 *                          ltfat_idgtreal_fb_plan_s** pout);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_init)(const LTFAT_REAL g[], const ltfatInt gl,
                             const ltfatInt a, const ltfatInt M, const dgt_phasetype ptype,
                             unsigned flags, LTFAT_NAME(idgtreal_fb_plan)** pout);

/** Execute plan for Inverse Discrete Gabor Transform for real signals using the filter bank algorithm
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
 * ltfat_idgtreal_fb_execute_d(ltfat_idgtreal_fb_plan_d* plan, const complex double c[],
 *                             const ltfatInt L, const ltfatInt W, double f[]);
 *
 * ltfat_idgtreal_fb_execute_s(ltfat_idgtreal_fb_plan_s* plan, const complex float c[],
 *                             const ltfatInt L, const ltfatInt W, float f[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_execute)(LTFAT_NAME(idgtreal_fb_plan)* p, const LTFAT_COMPLEX c[],
                                const ltfatInt L, const ltfatInt W, LTFAT_REAL f[]);

/** Destroy the plan
 *
 * \param[in]  plan   DGT plan
 * \return Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_fb_done_d(ltfat_idgtreal_fb_plan_d** plan);
 *
 * ltfat_idgtreal_fb_done_s(ltfat_idgtreal_fb_plan_s** plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_fb_done)(LTFAT_NAME(idgtreal_fb_plan)** p);

/** @}*/
