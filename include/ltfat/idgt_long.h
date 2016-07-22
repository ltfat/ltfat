typedef struct LTFAT_NAME(idgt_long_plan) LTFAT_NAME(idgt_long_plan);

/**
 *  \addtogroup dgtlong
 * @{
 *
 *
 */


/** Computes inverse Discrete Gabor Transform using the factorization algorithm
 *
 * \param[in]       c   Input coefficients, M x N x W array
 * \param[in]       g   Window to be used, array of length L
 * \param[in]       L   Signal length
 * \param[in]       W   Number of channels of the signal
 * \param[in]       a   Time hop factor
 * \param[in]       M   Number of frequency channels
 * \param[in]   ptype   Phase convention
 * \param[out]      f   Reconstructed signal, L x W array
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_long_d(const complex double c[], const double g[], const ltfatInt L,
 *                   const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                   const ltfat_phaseconvention ptype, complex double f[]);
 *
 * ltfat_idgt_long_s(const complex float c[], const float g[], const ltfatInt L,
 *                   const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                   const ltfat_phaseconvention ptype, complex float f[]);
 *
 * ltfat_idgt_long_dc(const complex double c[], const complex double g[], const ltfatInt L,
 *                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                    const ltfat_phaseconvention ptype, complex double f[]);
 *
 * ltfat_idgt_long_sc(const complex float c[], const complex float g[], const ltfatInt L,
 *                    const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                    const ltfat_phaseconvention ptype, complex float f[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgt_long)(const LTFAT_COMPLEX c[], const LTFAT_TYPE g[],
                      const ltfatInt L, const ltfatInt W,
                      const ltfatInt a, const ltfatInt M,
                      const ltfat_phaseconvention ptype, LTFAT_COMPLEX f[]);

/** Initialize inverse Discrete Gabor Transform plan for the factorization algorithm
 *
 * \note Please note that the input and output arrays will be overwritten when
 * anything else than FFTW_ESTIMATE is used in flags.
 *
 * \param[in]       c   Input coefficient array, M x N x W array
 * \param[in]       g   Window to be used, array of length L
 * \param[in]       L   Signal length
 * \param[in]       W   Number of channels of the signal
 * \param[in]       a   Time hop factor
 * \param[in]       M   Number of frequency channels
 * \param[in]       f   Reconstructed signal array, L x W array
 * \param[in]   flags   FFTW plan flags
 * \param[in]   ptype   Phase convention
 * \param[out]  pout    Initialized plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_long_init_d(complex double c[], const double g[], const ltfatInt L,
 *                        const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                        complex double f[], unsigned flags, const ltfat_phaseconvention ptype,
 *                        ltfat_idgt_long_plan_d** pout);
 *
 * ltfat_idgt_long_init_s(complex float c[], const float g[], const ltfatInt L,
 *                        const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                        complex float f[], unsigned flags, const ltfat_phaseconvention ptype,
 *                        ltfat_idgt_long_plan_s** pout);
 *
 * ltfat_idgt_long_init_dc(complex double c[], const complex double g[], const ltfatInt L,
 *                         const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                         complex double f[], unsigned flags, const ltfat_phaseconvention ptype,
 *                         ltfat_idgt_long_plan_dc** pout);
 *
 * ltfat_idgt_long_init_sc(complex float c[], const complex float g[], const ltfatInt L,
 *                         const ltfatInt W, const ltfatInt a, const ltfatInt M,
 *                         complex float f[], unsigned flags, const ltfat_phaseconvention ptype,
 *                         ltfat_idgt_long_plan_sc** pout);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgt_long_init)(LTFAT_COMPLEX c[], const LTFAT_TYPE g[],
                           const ltfatInt L, const ltfatInt W,
                           const ltfatInt a, const ltfatInt M,  LTFAT_COMPLEX f[],
                           unsigned flags, const ltfat_phaseconvention ptype, LTFAT_NAME(idgt_long_plan)** pout);


/** Execute the Inverse Discrete Gabor Transform plan
 *
 * \param[in]       p   Plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_long_execute_d(ltfat_idgt_long_plan_d* p);
 *
 * ltfat_idgt_long_execute_s(ltfat_idgt_long_plan_s* p);
 *
 * ltfat_idgt_long_execute_dc(ltfat_idgt_long_plan_dc* p);
 *
 * ltfat_idgt_long_execute_sc(ltfat_idgt_long_plan_sc* p);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgt_long_execute)(LTFAT_NAME(idgt_long_plan)* p);


/** Execute the Inverse Discrete Gabor Transform plan
 *
 * ... on arrays which might not have been used in init.
 *
 * \param[in]       p   Plan
 * \param[in]       c   Coefficients, size M x N xW
 * \param[out]      f   Output signal, size L x W  
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_long_execute_newarray_d(ltfat_idgt_long_plan_d* p, const complex double c[],
 *                                    complex double f[]);
 *
 * ltfat_idgt_long_execute_newarray_s(ltfat_idgt_long_plan_s* p, const complex double c[],
 *                                    complex double f[]);
 *
 * ltfat_idgt_long_execute_newarray_dc(ltfat_idgt_long_plan_dc* p, const complex double c[],
 *                                     complex double f[]);
 *
 * ltfat_idgt_long_execute_newarray_sc(ltfat_idgt_long_plan_sc* p, const complex double c[],
 *                                     complex double f[]);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgt_long_execute_newarray)(LTFAT_NAME(idgt_long_plan)* p,
                                       const LTFAT_COMPLEX c[],
                                       LTFAT_COMPLEX f[]);


/** Destroy the Discrete Gabor Transform plan
 *
 * \param[in]       plan   Plan
 * \returns Status code
 *
 * #### Versions #
 * <tt>
 * ltfat_idgt_long_done_d(ltfat_idgt_long_plan_d** p);
 *
 * ltfat_idgt_long_done_s(ltfat_idgt_long_plan_s** p);
 *
 * ltfat_idgt_long_done_dc(ltfat_idgt_long_plan_dc** p);
 *
 * ltfat_idgt_long_done_sc(ltfat_idgt_long_plan_sc** p);
 * </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(idgt_long_done)(LTFAT_NAME(idgt_long_plan)** plan);
/** @}*/

LTFAT_EXTERN void
LTFAT_NAME(idgt_fac)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *gf,
                     const ltfatInt L,
                     const ltfatInt W, const ltfatInt a, const ltfatInt M,
                     const ltfat_phaseconvention ptype, LTFAT_COMPLEX *f);

LTFAT_EXTERN void
LTFAT_NAME(idgt_walnut_execute)(LTFAT_NAME(idgt_long_plan)* p);
