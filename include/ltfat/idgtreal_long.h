typedef struct LTFAT_NAME(idgtreal_long_plan) LTFAT_NAME(idgtreal_long_plan);

/**
 *  \addtogroup dgtreallong
 * @{
 * For a detailed description see the dedicated page \ref dgttheory
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
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_d(const complex double c[], const double g[],
 *                       const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                       const ltfatInt M, const ltfat_phaseconvention ptype, double f[]);
 *
 * ltfat_idgtreal_long_s(const complex float c[], const float g[],
 *                       const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                       const ltfatInt M, const ltfat_phaseconvention ptype, float f[]);
 * </tt>
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | At least one of the following was NULL: \a f, \a g, \a c
 * LTFATERR_BADSIZE         | length of the signal and of the window \a L was less or equal to 0.
 * LTFATERR_NOTPOSARG       | At least one of the following was less or equal to zero: \a W, \a a, \a M
 * LTFATERR_BADTRALEN       | \a L is not divisible by both \a a and \a M.
 * LTFATERR_INITFAILED      | FFTW plan creation failed
 * LTFATERR_CANNOTHAPPEN    | \a ptype does not have a valid value from the ltfat_phaseconvention enum
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(idgtreal_long)(const LTFAT_COMPLEX c[], const LTFAT_REAL g[],
                          const ltfatInt L, const ltfatInt W,
                          const ltfatInt a, const ltfatInt M,
                          const ltfat_phaseconvention ptype, LTFAT_REAL f[]);

/** Initialize plan for Inverse Discrete Gabor Transform for real signals for the factorization algorithm
 *
 * \note \a f can be NULL if the plan is intended to be used with the _newarray execute function.
 * Similarly, \a c can also be NULL, but only if FFTW_ESTIMATE is passed in \a flags
 * (the FFTW planning routine does not touch the array). On the other hand,
 * the content of \a c might get overwritten if other FFTW planning flags are used.
 *
 * \param[in]     c   DGT coefficients, size M2 x N x W or NULL if (flags & FFTW_ESTIMATE) is nonzero.
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in]     f   Output signal, size L x W or NULL
 * \param[in] ptype   Phase convention
 * \param[in] flags   FFTW plan flags
 * \param[out] plan   IDGT plan
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_init_d(complex double c[], const double g[],
 *                            const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                            const ltfatInt M, double f[], const ltfat_phaseconvention ptype,
 *                            unsigned flags, ltfat_idgtreal_long_plan_d** plan);
 *
 * ltfat_idgtreal_long_init_s(complex float c[], const float g[],
 *                            const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                            const ltfatInt M, float f[], const ltfat_phaseconvention ptype,
 *                            unsigned flags, ltfat_idgtreal_long_plan_s** plan);
 * </tt>
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | The plan pointer \a plan or window \a g were NULL or \a c was NULL and (flags & FFTW_ESTIMATE) is zero
 * LTFATERR_BADSIZE         | Signal and window length \a L is less or equal to 0.
 * LTFATERR_NOTPOSARG       | Either of \a W, \a a, \a M was less or equal to zero.
 * LTFATERR_BADTRALEN       | \a L is not divisible by both \a a and \a M.
 * LTFATERR_INITFAILED      | The FFTW plan creation failed
 * LTFATERR_CANNOTHAPPEN    | \a ptype does not have a valid value from the ltfat_phaseconvention enum
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(idgtreal_long_init)(LTFAT_COMPLEX c[], const LTFAT_REAL g[],
                               const ltfatInt L, const ltfatInt W,
                               const ltfatInt a, const ltfatInt M, LTFAT_REAL f[],
                               const ltfat_phaseconvention ptype, unsigned flags,
                               LTFAT_NAME(idgtreal_long_plan)** plan);

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
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | The \a plan was NULL or it was created with \a f == NULL or \a c == NULL
 */
LTFAT_API int
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
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | Al least one of the arguments was NULL.
 */
LTFAT_API int
LTFAT_NAME(idgtreal_long_execute_newarray)(LTFAT_NAME(idgtreal_long_plan)* p,
        const LTFAT_COMPLEX* c, LTFAT_REAL* f);

/** Destroy the plan
 *
 * \param[in]  plan   IDGT plan
 *
 * #### Versions #
 * <tt>
 * ltfat_idgtreal_long_done_d(ltfat_idgtreal_long_plan_d** plan);
 *
 * ltfat_idgtreal_long_done_s(ltfat_idgtreal_long_plan_s** plan);
 * </tt>
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | plan or *plan was NULL.
 */
LTFAT_API int
LTFAT_NAME(idgtreal_long_done)(LTFAT_NAME(idgtreal_long_plan)** plan);

/** @}*/

LTFAT_API void
LTFAT_NAME(idgtreal_walnut_execute)(LTFAT_NAME(idgtreal_long_plan)* p);
