#ifndef _LTFAT_DGT_COMMON_H
#define _LTFAT_DGT_COMMON_H
#include "ltfat/basicmacros.h"
/** \defgroup dgtcommon Discrete Gabor Transform Common 
 * \addtogroup dgtcommon
 * @{
 */

/** Discrete Gabor Transform Phase convention
 *
 * There are two commonly used phase conventions. 
 *
 * Frequency invariant:
 * \anchor dgtfreqinv
 *  \f[
 *  c(m,n) 
 *   = \sum_{l=0}^{L-1}\! f(l)
 *   \overline{g(l-na)} \me^{-\mi 2\pi l m/M } \,
 *  \f]
 *  This is commonly used by mathematicans. Effectivelly, each
 *  frequency channel is demodulated to the baseband.
 *
 * Time invatiant:
 * \anchor dgttimeinv
 *  \f[
 *  c(m,n) 
 *   = \sum_{l=0}^{L-1}\! f(l)
 *   \overline{g(l-na)} \me^{-\mi 2\pi (l-na) m/M } \,
 *  \f]
 *  This is commonly used by engineers. This is equivalent to
 *  a filterbank.
 *
 */



typedef enum
{
    LTFAT_FREQINV = 0,
    LTFAT_TIMEINV = 1
} ltfat_phaseconvention;

/** @}*/

#ifdef __cplusplus
extern "C"
{
#endif

LTFAT_API int
ltfat_phaseconvention_is_valid(ltfat_phaseconvention in);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
