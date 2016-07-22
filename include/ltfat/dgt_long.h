typedef struct LTFAT_NAME(dgt_long_plan) LTFAT_NAME(dgt_long_plan);

/** \defgroup dgtlong Discrete Gabor Transform for long windows -- dgt_long 
 *  \addtogroup dgtlong
 * @{
 *

 *
 */

/** Computes DGT
 *
 * \param[in]      f  Multi-channel input signal, size L x W
 * \param[in]      g  Window, size L
 * \param[in]      L  Input signal length
 * \param[in]      W  Number of channels
 * \param[in]      a  Hop factor
 * \param[in]      M  Number of frequency channels (FFT size)
 * \param[in]  ptype  Number of frame diagonal samples
 * \param[out]     c  Output DGT coefficients, size M x N x W
 *
 * \returns Status code
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_dgt_long_d(const double f[], const double g[], const ltfatInt L,
 *                   const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                   const ltfat_phaseconvention ptype, complex double c[]);
 *
 *  ltfat_dgt_long_s(const float f[], const float g[], const ltfatInt L,
 *                   const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                   const ltfat_phaseconvention ptype, complex float c[]);
 *
 *  ltfat_dgt_long_dc(const complex double f[], const complex double g[],
 *                    const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                    const ltfatInt M, const ltfat_phaseconvention ptype,
 *                    complex double c[]);
 *
 *  ltfat_dgt_long_sc(const complex float f[], const complex float g[],
 *                    const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                    const ltfatInt M, const ltfat_phaseconvention ptype,
 *                    complex float c[]);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long)(const LTFAT_TYPE f[], const LTFAT_TYPE g[],
                     const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                     const ltfatInt M, const ltfat_phaseconvention ptype,
                     LTFAT_COMPLEX c[]);

/** Inicialization of the DGT plan
 *
 * \param[in]      f  Multi-channel input signal, size L x W
 * \param[in]      g  Window, size L
 * \param[in]      L  Input signal length
 * \param[in]      W  Number of channels
 * \param[in]      a  Hop factor
 * \param[in]      M  Number of frequency channels (FFT size)
 * \param[in]   cout  Output DGT coefficients, size M x N x W
 * \param[in]  ptype  Number of frame diagonal samples
 * \param[in]  flags  FFTW planning flag
 * \param[out]     p  DGT plan
 *
 * \returns Status code
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_dgt_long_init_d(const double f[], const double g[], const ltfatInt L,
 *                        const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                        complex double c[], const ltfat_phaseconvention ptype,
 *                        unsigned flags, dgt_long_plan_d** p);
 *
 *  ltfat_dgt_long_init_s(const float f[], const float g[], const ltfatInt L,
 *                        const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                        complex float c[], const ltfat_phaseconvention ptype,
 *                        unsigned flags, dgt_long_plan_s** p);
 *
 *  ltfat_dgt_long_init_dc(const complex double f[], const complex double g[],
 *                         const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                         const ltfatInt M, complex double c[],
 *                         const ltfat_phaseconvention ptype, unsigned flags,
 *                         dgt_long_plan_dc** p);
 *
 *  ltfat_dgt_long_init_sc(complex float f[], const complex float g[],
 *                         const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                         const ltfatInt M, complex float c[],
 *                         const ltfat_phaseconvention ptype, unsigned flags,
 *                         dgt_long_plan_sc** p);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long_init)(const LTFAT_TYPE f[], const LTFAT_TYPE g[], const ltfatInt L,
                          const ltfatInt W, const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX c[],
                          const ltfat_phaseconvention ptype, unsigned flags,
                          LTFAT_NAME(dgt_long_plan)** p);

/** Execute DGT plan
 *
 * \param[in]     plan  DGT plan
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_dgt_long_execute_d(ltfat_dgt_long_plan_d* plan);
 *
 *  ltfat_dgt_long_execute_s(ltfat_dgt_long_plan_s* plan);
 *
 *  ltfat_dgt_long_execute_dc(ltfat_dgt_long_plan_dc* plan);
 *
 *  ltfat_dgt_long_execute_sc(ltfat_dgt_long_plan_sc* plan);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long_execute)(LTFAT_NAME(dgt_long_plan)* plan);

/** Execute DGT plan
 *
 * ... on arrays which might not have been used in init.
 *
 * \param[in]     plan  DGT plan
 * \param[in]        f  Input signal, size L x W
 * \param[out]       c  Coefficients, size M x N x W
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_dgt_long_execute_newarray_d(ltfat_dgt_long_plan_d* plan, 
 *                                    const double f[], complex double c[]);
 *
 *  ltfat_dgt_long_execute_newarray_s(ltfat_dgt_long_plan_s* plan,
 *                                    const float f[], complex float c[]);
 *
 *  ltfat_dgt_long_execute_newarray_dc(ltfat_dgt_long_plan_cd* plan,
 *                                     const complex double f[], complex double c[]);
 *
 *  ltfat_dgt_long_execute_newarray_sc(ltfat_dgt_long_plan_cs* plan,
 *                                     const complex float f[], complex float c[]);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long_execute_newarray)(LTFAT_NAME(dgt_long_plan)* plan,
                                      const LTFAT_TYPE f[], LTFAT_COMPLEX c[]);


/** Destroy DGT plan
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_dgt_long_done_d(ltfat_dgt_long_plan_d** plan);
 *
 *  ltfat_dgt_long_done_s(ltfat_dgt_long_plan_s** plan);
 *
 *  ltfat_dgt_long_done_dc(ltfat_dgt_long_plan_dc** plan);
 *
 *  ltfat_dgt_long_done_sc(ltfat_dgt_long_plan_sc** plan);
 *  </tt>
 */

LTFAT_EXTERN int
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan)** plan);

/** @}*/

LTFAT_EXTERN int
LTFAT_NAME(dgt_walnut_execute)(LTFAT_NAME(dgt_long_plan)* plan, LTFAT_COMPLEX* cout);
