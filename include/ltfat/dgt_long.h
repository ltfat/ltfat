typedef struct LTFAT_NAME(dgt_long_plan) LTFAT_NAME(dgt_long_plan);

/** \defgroup dgtlong Discrete Gabor Transform for long windows -- dgt_long 
 *  \addtogroup dgtlong
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
 *  dgt_long_d(const double *f, const double *g, const ltfatInt L,
 *             const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *             const dgt_phasetype ptype, complex double *c);
 *
 *  dgt_long_s(const float *f, const float *g, const ltfatInt L,
 *             const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *             const dgt_phasetype ptype, complex float *c);
 *
 *  dgt_long_cd(const complex double *f, const complex double *g,
 *              const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *              const ltfatInt M, const dgt_phasetype ptype,
 *              complex double *c);
 *
 *  dgt_long_cs(const complex float *f, const complex float *g,
 *              const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *              const ltfatInt M, const dgt_phasetype ptype,
 *              complex float *c);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                     const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                     const ltfatInt M, const dgt_phasetype ptype,
                     LTFAT_COMPLEX *c);

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
 *  dgt_long_init_d(const double *f, const double *g, const ltfatInt L,
 *                  const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                  complex double *c, const dgt_phasetype ptype,
 *                  unsigned flags, dgt_long_plan_d** p);
 *
 *  dgt_long_init_s(const float *f, const float *g, const ltfatInt L,
 *                  const ltfatInt W,  const ltfatInt a, const ltfatInt M,
 *                  complex float *c, const dgt_phasetype ptype,
 *                  unsigned flags, dgt_long_plan_s** p);
 *
 *  dgt_long_init_cd(const complex double *f, const complex double *g,
 *                   const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                   const ltfatInt M, complex double *c,
 *                   const dgt_phasetype ptype, unsigned flags,
 *                   dgt_long_plan_cd** p);
 *
 *  dgt_long_init_cs(complex float *f, const complex float *g,
 *                   const ltfatInt L, const ltfatInt W,  const ltfatInt a,
 *                   const ltfatInt M, complex float *c,
 *                   const dgt_phasetype ptype, unsigned flags,
 *                   dgt_long_plan_cs** p);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long_init)(const LTFAT_TYPE *f, const LTFAT_TYPE *g, const ltfatInt L,
                          const ltfatInt W, const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *c,
                          const dgt_phasetype ptype, unsigned flags,
                          LTFAT_NAME(dgt_long_plan)** p);

/** Execute DGT plan
 *
 * \param[in]     plan  DGT plan
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  dgt_long_execute_d(dgt_long_plan_d* plan);
 *
 *  dgt_long_execute_s(dgt_long_plan_s* plan);
 *
 *  dgt_long_execute_cd(dgt_long_plan_cd* plan);
 *
 *  dgt_long_execute_cs(dgt_long_plan_cs* plan);
 *  </tt>
 */
LTFAT_EXTERN int
LTFAT_NAME(dgt_long_execute)(LTFAT_NAME(dgt_long_plan)* plan);

/** Destroy DGT plan
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  dgt_long_done_d(dgt_long_plan_d** plan);
 *
 *  dgt_long_done_s(dgt_long_plan_s** plan);
 *
 *  dgt_long_done_cd(dgt_long_plan_cd** plan);
 *
 *  dgt_long_done_cs(dgt_long_plan_cs** plan);
 *  </tt>
 */

LTFAT_EXTERN int
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan)** plan);

/** @}*/

LTFAT_EXTERN int
LTFAT_NAME(dgt_walnut_execute)(LTFAT_NAME(dgt_long_plan)* plan);
