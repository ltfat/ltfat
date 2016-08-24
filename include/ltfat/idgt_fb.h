typedef struct LTFAT_NAME(idgt_fb_plan) LTFAT_NAME(idgt_fb_plan);

/** 
 *  \addtogroup dgtfb
 * @{
 * For a detailed description see the dedicated page \ref dgttheory
 */

/** Compute Inverse Discrete Gabor Transform using filter bank algorithm
 *
 * \param[out]    c   DGT coefficients, size M x N x W
 * \param[in]     g   Window, size L x 1
 * \param[in]     L   Signal length
 * \param[in]    gl   Window length
 * \param[in]     W   Number of channels of the signal
 * \param[in]     a   Time hop facto
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[in]     f   Output signal, size L x W
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_fb_d(const complex double c[], const double g[],
 *                 const ltfatInt L, const ltfatInt gl,
 *                 const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                 const ltfat_phaseconvention ptype, complex double f[]);
 *
 * ltfat_idgt_fb_s(const complex float c[], const float g[],
 *                 const ltfatInt L, const ltfatInt gl,
 *                 const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                 const ltfat_phaseconvention ptype, complex float f[]);
 *
 * ltfat_idgt_fb_dc(const complex double c[], const complex double g[],
 *                  const ltfatInt L, const ltfatInt gl,
 *                  const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                  const ltfat_phaseconvention ptype, complex double f[]);
 *
 * ltfat_idgt_fb_sc(const complex float c[], const complex float g[],
 *                  const ltfatInt L, const ltfatInt gl,
 *                  const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                  const ltfat_phaseconvention ptype, complex float f[]);
 * </tt>
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | At least one of the following was NULL: \a f, \a g, \a c
 * LTFATERR_BADSIZE         | Length of the signal \a L or the length of the window \a gl was less or equal to 0.
 * LTFATERR_NOTPOSARG       | At least one of the following was less or equal to zero: \a W, \a a, \a M
 * LTFATERR_BADTRALEN       | \a L must be bigger or equal to \a gl and must be divisible by \a a
 * LTFATERR_INITFAILED      | FFTW plan creation failed
 * LTFATERR_CANNOTHAPPEN    | \a ptype does not have a valid value from the ltfat_phaseconvention enum
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(idgt_fb)(const LTFAT_COMPLEX c[], const LTFAT_TYPE g[],
                    const ltfatInt L, const ltfatInt gl,
                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
                    const ltfat_phaseconvention ptype, LTFAT_COMPLEX f[]);

/** Initialize plan for Inverse Discrete Gabor Transform for the filter bank algorithm
 *
 * \param[in]     g   Window, size L x 1
 * \param[in]    gl   Window length
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in] ptype   Phase convention
 * \param[in] flags   FFTW plan flags
 * \param[out] plan   DGT plan
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_fb_init_d(const double g[], const ltfatInt gl, const ltfatInt a,
 *                      const ltfatInt M, const ltfat_phaseconvention ptype, unsigned flags
 *                      ltfat_idgt_fb_plan_d** plan);
 *
 * ltfat_idgt_fb_init_s(const float g[], const ltfatInt gl, const ltfatInt a,
 *                      const ltfatInt M, const ltfat_phaseconvention ptype, unsigned flags
 *                      ltfat_idgt_fb_plan_s** plan);
 *
 * ltfat_idgt_fb_init_dc(const complex double g[], const ltfatInt gl, const ltfatInt a,
 *                       const ltfatInt M, const ltfat_phaseconvention ptype, unsigned flags
 *                       ltfat_idgt_fb_plan_dc** plan);
 *
 * ltfat_idgt_fb_init_sc(const complex float g[], const ltfatInt gl, const ltfatInt a,
 *                       const ltfatInt M, const ltfat_phaseconvention ptype, unsigned flags
 *                       ltfat_idgt_fb_plan_sc** plan);
 * </tt>
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | At least one of the following was NULL: \a g, \a plan
 * LTFATERR_BADSIZE         | Length of the window \a gl was less or equal to 0.
 * LTFATERR_NOTPOSARG       | At least one of the following was less or equal to zero: \a a, \a M
 * LTFATERR_INITFAILED      | FFTW plan creation failed
 * LTFATERR_CANNOTHAPPEN    | \a ptype does not have a valid value from the ltfat_phaseconvention enum
 * LTFATERR_NOMEM           | Indicates that heap allocation failed
 */
LTFAT_API int
LTFAT_NAME(idgt_fb_init)(const LTFAT_TYPE g[], const ltfatInt gl,
                         const ltfatInt a, const ltfatInt M, const ltfat_phaseconvention ptype,
                         unsigned flags, LTFAT_NAME(idgt_fb_plan)** plan);

/** Execute plan for Inverse Discrete Gabor Transform using the filter bank algorithm
 *
 * \param[in]  plan   DGT plan
 * \param[in]     f   Input signal, size L x W
 * \param[in]     L   Signal length
 * \param[in]     W   Number of channels of the signal
 * \param[out]    c   DGT coefficients, size M x N x W
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_fb_execute_d(ltfat_idgt_fb_plan_d* plan, const complex double c[],
 *                         const ltfatInt L, const ltfatInt W, double f[]);
 *
 * ltfat_idgt_fb_execute_s(ltfat_idgt_fb_plan_s* plan, const complex float c[],
 *                         const ltfatInt L, const ltfatInt W, float f[]);
 *
 * ltfat_idgt_fb_execute_dc(ltfat_idgt_fb_plan_dc* plan, const complex double c[],
 *                          const ltfatInt L, const ltfatInt W, double f[]);
 *
 * ltfat_idgt_fb_execute_sc(ltfat_idgt_fb_plan_sc* plan, const complex float c[],
 *                          const ltfatInt L, const ltfatInt W, float f[]);
 * </tt>
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | At least one of the following was NULL: \a f, \a c, \a plan
 * LTFATERR_BADSIZE         | Length of the signal \a L was less or equal to 0.
 * LTFATERR_BADTRALEN       | \a L must be bigger or equal to \a gl and must be divisible by \a a
 * LTFATERR_NOTPOSARG       | \a W was less or equal to 0.
 */
LTFAT_API int
LTFAT_NAME(idgt_fb_execute)(LTFAT_NAME(idgt_fb_plan)* p, const LTFAT_COMPLEX c[],
                            const ltfatInt L, const ltfatInt W, LTFAT_COMPLEX f[]);

/** Destroy the plan
 *
 * \param[in]  plan   DGT plan
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_fb_done_d(ltfat_idgt_fb_plan_d** plan);
 *
 * ltfat_idgt_fb_done_s(ltfat_idgt_fb_plan_s** plan);
 *
 * ltfat_idgt_fb_done_dc(ltfat_idgt_fb_plan_dc** plan);
 *
 * ltfat_idgt_fb_done_sc(ltfat_idgt_fb_plan_sc** plan);
 * </tt>
 *
 * \returns
 * Status code              | Description
 * -------------------------|--------------------------------------------
 * LTFATERR_SUCCESS         | Indicates no error
 * LTFATERR_NULLPOINTER     | plan or *plan was NULL.
 */
LTFAT_API int
LTFAT_NAME(idgt_fb_done)(LTFAT_NAME(idgt_fb_plan)** p);


/** @}*/
