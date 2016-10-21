/** \defgroup windows Gabor Windows
 * \addtogroup windows
 * @{
 */

#ifndef _LTFAT_CI_WINDOWS_H
#define _LTFAT_CI_WINDOWS_H

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
    LTFAT_TRUNCGAUSS01
}
LTFAT_FIRWIN;

#endif /* _CI_WINDOWS_H */



/** Creates real, whole-point symmetric, zero delay window.
 *
 * \note Please note that the window format is slightly unusual i.e.
 * the (unique) peak of the window is at index [0] of the array and the 
 * left tail is circularly wrapped such that the ends of both tails 
 * meet in the middle of the array. 
 * fftshift() can be used to format the array to the usual format
 * with peak in the middle. 
 *
 * \see normalize fftshift
 *
 * \param[in]   win  Window type
 * \param[in]   gl   Window length
 * \param[out]  g    Window
 *
 * #### Function versions #
 * <tt>
 * ltfat_firwin_d(LTFAT_FIRWIN win, ltfat_int gl, double* g);
 *
 * ltfat_firwin_s(LTFAT_FIRWIN win, ltfat_int gl, float* g);
 *
 * ltfat_firwin_dc(LTFAT_FIRWIN win, ltfat_int gl, ltfat_complex_d* g);
 *
 * ltfat_firwin_sc(LTFAT_FIRWIN win, ltfat_int gl, ltfat_complex_s* g);
 * </tt>
 * \returns
 * Status code           | Description
 * ----------------------|--------------------------------------------
 * LTFATERR_SUCCESS      | Indicates no error
 * LTFATERR_NULLPOINTER  | The output array is NULL
 * LTFATERR_BADSIZE      | Length of the array is less or equal to 0.
 * LTFATERR_CANNOTHAPPEN | \a win is not a valid value from the LTFAT_FIRWIN enum
 */
LTFAT_API int
LTFAT_NAME(firwin)(LTFAT_FIRWIN win, ltfat_int gl, LTFAT_TYPE* g);


/** @} */
