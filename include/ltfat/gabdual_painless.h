/** \defgroup gabdual Canonical dual and tight windows
 */

/** Compute the first dl samples of the Gabor frame operator diagonal
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[in]  dl    Number of frame diagonal samples
 * \param[out]  d    Frame diagonal
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(gabframediag)(const LTFAT_TYPE* g, ltfatInt gl,
                         ltfatInt a, ltfatInt M, ltfatInt dl, LTFAT_TYPE* d);

/** \addtogroup gabdual
 * @{
 */

/** Compute canonical dual window for painless Gabor system
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gd    Canonical dual window
 *
 * #### Versions #
 * <tt>
 * ltfat_gabdual_painless_d(const double g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M, double gd[]);
 *
 * ltfat_gabdual_painless_s(const float g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M, float gd[]);
 *
 * ltfat_gabdual_painless_cd(const complex double g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M,
 *                           complex double gd[]);
 *
 * ltfat_gabdual_painless_cs(const complex float g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M,
 *                           complex float gd[]);
 * </tt>
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL.
 * LTFATERR_BADSIZE      | Length of the windows is less or equal to 0.
 * LTFATERR_NOTPOSARG    | Either of \a a, \a M is less or equal to zero.
 * LTFATERR_NOTAFRAME    | System is not a frame.
 * LTFATERR_NOTPAINLESS  | System is not painless. 
 * LTFATERR_NOMEM        | Indicates that the heap allocation failed.
 */
LTFAT_EXTERN int
LTFAT_NAME(gabdual_painless)(const LTFAT_TYPE g[], const ltfatInt gl, const ltfatInt a,
                             const ltfatInt M, LTFAT_TYPE gd[]);

/** Compute canonical tight window for painless Gabor system
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gt    Canonical tight window
 *
 * #### Versions #
 * <tt>
 * ltfat_gabtight_painless_d(const double g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M, double gt[]);
 *
 * ltfat_gabtight_painless_s(const float g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M, float gt[]);
 *
 * ltfat_gabtight_painless_cd(const complex double g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M,
 *                            complex double gt[]);
 *
 * ltfat_gabtight_painless_cs(const complex float g[], const ltfatInt gl, const ltfatInt a, const ltfatInt M,
 *                            complex float gt[]);
 * </tt>
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL.
 * LTFATERR_BADSIZE      | Length of the windows is less or equal to 0.
 * LTFATERR_NOTPOSARG    | Either of \a a, \a M is less or equal to zero.
 * LTFATERR_NOTAFRAME    | System is not a frame.
 * LTFATERR_NOTPAINLESS  | System is not painless. 
 * LTFATERR_NOMEM        | Indicates that the heap allocation failed.
 */
LTFAT_EXTERN int
LTFAT_NAME(gabtight_painless)(const LTFAT_TYPE g[], const ltfatInt gl, const ltfatInt a,
                              const ltfatInt M, LTFAT_TYPE gt[]);


/** @} */

/** Compute partition of unity window from any window
 *
 * Creates a window which is PU for given lattice.
 *
 * When using PU window for either analysis or synthesis, the other window can just be
 * rectangular.
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gpu   Partition of unity window
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(gabpu_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                           ltfatInt M, LTFAT_TYPE* gpu);
