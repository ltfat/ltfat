
/** \addtogroup utils
 * @{
 */

/** Circshift in the Fourier domain via modulation
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[in]  shift  Shift amount (can be non-integer)
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_fftcircshift_d(const complex double in[], const ltfatInt L,const double shift, complex double out[]);
 *
 *  ltfat_fftcircshift_s(const complex float in[], const ltfatInt L,const double shift, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(fftcircshift)(const LTFAT_COMPLEX in[], const ltfatInt L, const double shift,
                         LTFAT_COMPLEX out[]);

/** fftshift in the Fourier domain via modulation
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_fftfftshift_d(const complex double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_fftfftshift_s(const complex float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(fftfftshift)(const LTFAT_COMPLEX in[], const ltfatInt L, LTFAT_COMPLEX out[]);

/** ifftshift in the Fourier domain via modulation
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_fftifftshift_d(const complex double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_fftifftshift_s(const complex float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(fftifftshift)(const LTFAT_COMPLEX in[], const ltfatInt L, LTFAT_COMPLEX out[]);

/** circshift in the Fourier domain (fftreal) via modulation
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[in]  shift  Shift amount (can be non-integer)
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_fftrealcircshift_d(const complex double in[], const ltfatInt L,const double shift, complex double out[]);
 *
 *  ltfat_fftrealcircshift_s(const complex float in[], const ltfatInt L,const double shift, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(fftrealcircshift)( const LTFAT_COMPLEX in[], const ltfatInt L, const double shift,
                              LTFAT_COMPLEX out[]);

/** fftshift in the Fourier domain (fftreal) via modulation
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_fftrealfftshift_d(const complex double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_fftrealfftshift_s(const complex float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(fftrealfftshift)(const LTFAT_COMPLEX in[], const ltfatInt L, LTFAT_COMPLEX out[]);

/** ifftshift in the Fourier domain (fftreal) via modulation
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_fftrealifftshift_d(const complex double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_fftrealifftshift_s(const complex float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */

LTFAT_EXTERN int
LTFAT_NAME(fftrealifftshift)(const LTFAT_COMPLEX in[], const ltfatInt L, LTFAT_COMPLEX out[]);


/** Real to complex array
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltfat_real2complex_array_d(const double in[], const ltfatInt L, complex double out[]);
 *
 *  ltfat_real2complex_array_s(const float in[], const ltfatInt L, complex float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(real2complex_array)(const LTFAT_REAL in[], const ltfatInt L, LTFAT_COMPLEX out[]);

/** Complex to real array
 *
 * \param[in]     in  Input array
 * \param[in]      L  Length of arrays
 * \param[out]   out  Output array
 *
 *  Function versions
 *  -----------------
 *
 *  <tt>
 *  ltaft_complex2real_array_d(const complex double in[], const ltfatInt L, double out[]);
 *
 *  ltfat_complex2real_array_s(const complex float in[], const ltfatInt L, float out[]);
 *  </tt>
 *
 * \returns Status code
 */
LTFAT_EXTERN int
LTFAT_NAME(complex2real_array)(const LTFAT_COMPLEX in[], const ltfatInt L, LTFAT_REAL out[]);


/** @}*/
