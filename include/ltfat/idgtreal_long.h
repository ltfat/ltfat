typedef struct LTFAT_NAME(idgtreal_long_plan) LTFAT_NAME(idgtreal_long_plan);

/**
 *  \addtogroup dgtreallong
 * @{
 *
 *
 */

/** Compute Inverse Discrete Gabor Transform for real signals using the factorization algorithm
 *
 * \param[in]     c   DGT coefficients, size M2 x N x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[out]    f   Output signal, size L x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_d(const complex double c[], const double g[],
 *                       const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                       const ltfatInt M, const dgt_phasetype ptype, double f[]);
 *
 * ltfat_idgtreal_long_s(const complex float c[], const float g[],
 *                       const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                       const ltfatInt M, const dgt_phasetype ptype, float f[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long)(const LTFAT_COMPLEX c[], const LTFAT_REAL g[],
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M,
                          const dgt_phasetype ptype, LTFAT_REAL f[]);

/** Initialize plan for Inverse Discrete Gabor Transform for real signals for the factorization algorithm
 *
 * \param[in]     c   DGT coefficients, size M2 x N x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in]     f   Output signal, size L x W
 * \param[in] ptype   Phase convention
 * \param[in] flags   FFTW plan flags
 * \param[out] pout   IDGT plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_init_d(complex double c[], const double g[],
 *                            const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                            const ltfatInt M, double f[], const dgt_phasetype ptype,
 *                            unsigned flags, ltfat_idgtreal_long_plan_d** pout);
 *
 * ltfat_idgtreal_long_init_s(complex float c[], const float g[],
 *                            const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                            const ltfatInt M, float f[], const dgt_phasetype ptype,
 *                            unsigned flags, ltfat_idgtreal_long_plan_s** pout);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_init)(LTFAT_COMPLEX c[], const LTFAT_REAL g[],
                               const ltfatInt L, const ltfatInt W,
                               const ltfatInt a, const ltfatInt M, LTFAT_REAL f[],
                               const dgt_phasetype ptype, unsigned flags,
                               LTFAT_NAME(idgtreal_long_plan)** pout);

/** Execute plan for Inverse Discrete Gabor Transform for real signals using the factorization algorithm
 *
 * \param[in]  plan   IDGT plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_execute_d(ltfat_idgtreal_long_plan_d* plan);
 *
 * ltfat_idgtreal_long_execute_s(ltfat_idgtreal_long_plan_s* plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_execute)(LTFAT_NAME(idgtreal_long_plan)* plan);

/** Execute plan for Inverse Discrete Gabor Transform for real signals using the factorization algorithm
 *
 * ... on arrays which might not have been used in init.
 *
 * \param[in]  plan   IDGT plan
 * \param[in]     c   Coefficients, size M2 x N x W  
 * \param[out]    f   Output signal, size L x W
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_execute_newarray_d(ltfat_idgtreal_long_plan_d* plan,
 *                                        const complex double c[], double f[]);
 *
 * ltfat_idgtreal_long_execute_newarray_s(ltfat_idgtreal_long_plan_s* plan,
 *                                        const complex float c[], float f[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_execute_newarray)(LTFAT_NAME(idgtreal_long_plan)* p,
        const LTFAT_COMPLEX* c, LTFAT_REAL* f);

/** Destroy the plan
 *
 * \param[in]  plan   IDGT plan
 * \return Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_done_d(ltfat_idgtreal_long_plan_d** plan);
 *
 * ltfat_idgtreal_long_done_s(ltfat_idgtreal_long_plan_s** plan);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgtreal_long_done)(LTFAT_NAME(idgtreal_long_plan)** plan);

/** @}*/

LTFAT_EXTERN void
LTFAT_NAME(idgtreal_walnut_execute)(LTFAT_NAME(idgtreal_long_plan)* p);
