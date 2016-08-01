/** \addtogroup utils
 * @{
 * Utility functions
 *
 * \note All utility functions working with arrays can work inplace (i.e. in == out) without
 * allocating additional memory internally.
 */

#ifndef _ciutils_h
#define _ciutils_h

typedef enum
{
    /** Don't normalize */
    LTFAT_NORMALIZE_NULL = 0,
    /** Normalize to 1 norm (divide by the sum of abs. values)  */
    /**@{*/
    LTFAT_NORMALIZE_AREA,
    LTFAT_NORMALIZE_1 = LTFAT_NORMALIZE_AREA,
    /**@}*/
    /** Normalize to 2 norm (divide by the square root of sum of squares of abs. values) */
    /**@{*/
    LTFAT_NORMALIZE_ENERGY,
    LTFAT_NORMALIZE_2 = LTFAT_NORMALIZE_ENERGY,
    /**@}*/
    /** Normalize to inf norm (divide by the max abs. val.)*/
    /**@{*/
    LTFAT_NORMALIZE_INF,
    LTFAT_NORMALIZE_PEAK = LTFAT_NORMALIZE_INF,
    /**@}*/
    // LTFAT_NORMALIZE_RMS,
    // LTFAT_NORMALIZE_WAV,
} ltfat_normalize_t;

#endif

/** Shift array circurarly
 *
 *  Works exactly like
 *  <a href="http://de.mathworks.com/help/matlab/ref/circshift.html">circshift</a>
 *  form Matlab.
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[in]  shift  Shift amount (can be negative)
 * \param[out]   out  Output array
 *
 *  #### Function versions ####
 *  <tt>
 *  ltfat_circshift_d(const double in[], const ltfatInt L,const ltfatInt shift, double out[]);
 *
 *  ltfat_circshift_s(const float in[], const ltfatInt L,const ltfatInt shift, float out[]);
 *
 *  ltfat_circshift_dc(const complex double in[], const ltfatInt L,const ltfatInt shift, complex double out[]);
 *
 *  ltfat_circshift_sc(const complex float in[], const ltfatInt L,const ltfatInt shift, complex float out[]);
 *  </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 */
LTFAT_EXTERN int
LTFAT_NAME(circshift)(const LTFAT_TYPE in[], const ltfatInt L,
                      const ltfatInt shift, LTFAT_TYPE out[]);

/** fftshift an array
 *
 *  Works exactly like
 *  <a href="http://de.mathworks.com/help/matlab/ref/fftshift.html">fftshift</a>
 *  form Matlab for vectors.
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  #### Function versions ####
 *  <tt>
 *  ltfat_fftshift_d(const double in[], const ltfatInt L, double out[]);
 *
 *  ltfat_fftshift_s(const float in[], const ltfatInt L, float out[]);
 *
 *  ltfat_fftshift_dc(const complex double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_fftshift_sc(const complex float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 */
LTFAT_EXTERN int
LTFAT_NAME(fftshift)(const LTFAT_TYPE in[], ltfatInt L, LTFAT_TYPE out[]);

/** ifftshift an array
 *
 *  Works exactly like
 *  <a href="http://de.mathworks.com/help/matlab/ref/ifftshift.html">ifftshift</a>
 *  form Matlab for vectors. Undoes the action of fftshift
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  #### Function versions ####
 *  <tt>
 *  ltfat_ifftshift_d(const double in[], const ltfatInt L, double out[]);
 *
 *  ltfat_ifftshift_s(const float in[], const ltfatInt L, float out[]);
 *
 *  ltfat_ifftshift_dc(const complex double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_ifftshift_sc(const complex float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 */
LTFAT_EXTERN int
LTFAT_NAME(ifftshift)(const LTFAT_TYPE in[], ltfatInt L, LTFAT_TYPE out[]);


/** Extend FIR window to long window
 *
 *  Works exactly like
 *  <a href="http://ltfat.github.io/doc/sigproc/fir2long.html">fir2long</a>
 *  form LTFAT i.e. extends \a in by inserting zeros in the middle.
 *  \a Llong must be greater or equal to \a Lfir.
 *
 * \param[in]     in  Input array
 * \param[in]   Lfir  Length of input array
 * \param[in]  Llong  Length of output array
 * \param[out]   out  Output array
 *
 *  #### Function versions ####
 *  <tt>
 *  ltfat_fir2long_d(const double in[], const ltfatInt Lfir, const ltfatInt Llong, double out[]);
 *
 *  ltfat_fir2long_s(const float in[], const ltfatInt Lfir, const ltfatInt Llong, float out[]);
 *
 *  ltfat_fir2long_dc(const complex double in[], const ltfatInt Lfir, const ltfatInt Llong, complex double out[]);
 *
 *  ltfat_fir2long_sc(const complex float in[], const ltfatInt Lfir, const ltfatInt Llong, complex float out[]);
 *  </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 * LTFATERR_BADREQSIZE   | Output array is shorter than the input array: \a Llong < \a Lfir
 */
LTFAT_EXTERN int
LTFAT_NAME(fir2long)(const LTFAT_TYPE in[], const ltfatInt Lfir, const ltfatInt Llong,
                     LTFAT_TYPE out[]);

/** Cut long window to a FIR window
 *
 *  Works exactly like
 *  <a href="http://ltfat.github.io/doc/sigproc/long2fir.html">long2fir</a>
 *  form LTFAT i.e. it removes the middle part.
 *  Llong must be greater or equal to Lfir.
 *
 * \param[in]     in  Input array
 * \param[in]  Llong  Length of input array
 * \param[in]   Lfir  Length of output array
 * \param[out]   out  Output array
 *
 *  #### Function versions ####
 *  <tt>
 *  ltfat_long2fir_d(const double in[], const ltfatInt Llong, const ltfatInt Lfir, double out[]);
 *
 *  ltfat_long2fir_s(const float in[], const ltfatInt Llong, const ltfatInt Lfir, float out[]);
 *
 *  ltfat_long2fir_dc(const complex double in[], const ltfatInt Llong, const ltfatInt Lfir, complex double out[]);
 *
 *  ltfat_long2fir_sc(const complex float in[], const ltfatInt Llong, const ltfatInt Lfir, complex float out[]);
 *  </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 * LTFATERR_BADREQSIZE   | Output array is longer than the input array: \a Lfir > \a Llong
 */
LTFAT_EXTERN int
LTFAT_NAME(long2fir)(const LTFAT_TYPE in[], const ltfatInt Llong, const ltfatInt Lfir,
                     LTFAT_TYPE out[]);

/** Normalize a vector
 *
 * Normalizes the input array such that the chosen norm is 1.
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of the arrays
 * \param[in]   flag  Norm
 * \param[out]   out  Output array
 *
 *  #### Function versions ####
 *  <tt>
 *  ltfat_normalize_d(const double in[], const ltfatInt L, ltfat_normalize_t flag, double out[]);
 *
 *  ltfat_normalize_s(const float in[], const ltfatInt L, ltfat_normalize_t flag, float out[]);
 *
 *  ltfat_normalize_dc(const complex double in[], const ltfatInt L, ltfat_normalize_t flag, complex double out[]);
 *
 *  ltfat_normalize_sc(const complex float in[], const ltfatInt L, ltfat_normalize_t flag, complex float out[]);
 *  </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 * LTFATERR_CANNOTHAPPEN | \a flag is not defined in ltfat_normalize_t enum.
 */
LTFAT_EXTERN int
LTFAT_NAME(normalize)(const LTFAT_TYPE in[], const ltfatInt L,
                      ltfat_normalize_t flag, LTFAT_TYPE out[]);

/** Ensure the array has complex interleaved layout
 *
 * This is a convenience function.
 * Obviously the *_dc and *_sc versions of the function do nothing.
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 * #### Function versions ####
 * <tt>
 * ltfat_ensurecomplex_array_d(const double in[], const ltfatInt L, complex double out[]);
 *
 * ltfat_ensurecomplex_array_s(const float in[], const ltfatInt L, complex float out[]);
 *
 * ltfat_ensurecomplex_array_dc(const complex double in[], const ltfatInt L, complex double out[]);
 *
 * ltfat_ensurecomplex_array_sc(const complex float in[], const ltfatInt L, complex float out[]);
 * </tt>
 *
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | Either of the arrays is NULL
 * LTFATERR_BADSIZE      | Length of the arrays is less or equal to 0.
 */
LTFAT_EXTERN int
LTFAT_NAME(ensurecomplex_array)(const LTFAT_TYPE *in, const ltfatInt L, LTFAT_COMPLEX *out);

/** @}*/

LTFAT_EXTERN void
LTFAT_NAME(dgtphaselockhelper)(LTFAT_TYPE *cin, const ltfatInt L,
                               const ltfatInt W, const ltfatInt a,
                               const ltfatInt M, LTFAT_TYPE *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgtphaseunlockhelper)(LTFAT_TYPE *cin, const ltfatInt L,
                                 const ltfatInt W, const ltfatInt a,
                                 const ltfatInt M, LTFAT_TYPE *cout);

LTFAT_EXTERN int
LTFAT_NAME(reverse_array)(const LTFAT_TYPE *in, const ltfatInt L, LTFAT_TYPE *out);

LTFAT_EXTERN int
LTFAT_NAME(conjugate_array)(const LTFAT_TYPE *in, const ltfatInt L, LTFAT_TYPE *out);

LTFAT_EXTERN int
LTFAT_NAME(periodize_array)(const LTFAT_TYPE *in, const ltfatInt Lin,
                            const ltfatInt Lout, LTFAT_TYPE *out );

LTFAT_EXTERN int
LTFAT_NAME(fold_array)(const LTFAT_TYPE *in, const ltfatInt Lin,
                       const ltfatInt offset,
                       const ltfatInt Lfold, LTFAT_TYPE *out);

LTFAT_EXTERN
void LTFAT_NAME(findmaxinarray)(const LTFAT_TYPE *in, const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx);

LTFAT_EXTERN int
LTFAT_NAME(findmaxinarraywrtmask)(const LTFAT_TYPE *in, const int *mask, const ltfatInt L, LTFAT_TYPE* max, ltfatInt* idx);
