/** \addtogroup dgt
 *  @{
 *
 * @file gabdualfir.h
 * @author Zdeněk Průša
 * @date 2 June 2016
 * @brief GABDUALFIR header
 *
 */

#ifndef _GABDUALFIR_H
#define _GABDUALFIR_H

typedef enum
{
    LTFAT_HANN, LTFAT_HANNING, LTFAT_NUTTALL10,
    LTFAT_SQRTHANN, LTFAT_COSINE, LTFAT_SINE,
    LTFAT_HAMMING,
    LTFAT_NUTTALL01,
    LTFAT_SQUARE, LTFAT_RECT,
    LTFAT_TRIA, LTFAT_TRIANGULAR, LTFAT_BARTLETT,
    LTFAT_SQRTTRIA,
    LTFAT_BLACKMAN,
    LTFAT_BLACKMAN2,
    LTFAT_NUTTALL, LTFAT_NUTTALL12,
    LTFAT_OGG, LTFAT_ITERSINE,
    LTFAT_NUTTALL20,
    LTFAT_NUTTALL11,
    LTFAT_NUTTALL02,
    LTFAT_NUTTALL30,
    LTFAT_NUTTALL21,
    LTFAT_NUTTALL03,
}
LTFAT_FIRWIN;

#endif /* _GABDUALFIR_H */

LTFAT_EXTERN int
LTFAT_NAME(firwin)(LTFAT_FIRWIN win, int gl, LTFAT_TYPE* g);

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

/** Compute canonical tight window for painless Gabor system
 *
 * Painless system requires gl<=M and gl>=a
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gt    Canonical tight window
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(gabtight_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                              ltfatInt M, LTFAT_TYPE* gt);

/** Compute canonical dual window for painless Gabor system
 *
 * Painless system requires gl<=M and gl>=a
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gd    Canonical dual window
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(gabdual_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a, ltfatInt M,
                             LTFAT_TYPE* gd);

/** Compute partition of unity window from any window
 *
 * Creates a window which is PU for given lattice.
 * 
 * When using PU window for either analysis or synthesis, the other window can just be
 * rectangular.
 *
 * Painless system requires gl<=M and gl>=a
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

