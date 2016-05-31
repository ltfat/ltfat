/** \addtogroup rtdgtreal
 * @{
 *  RTDGTREAL - Real-Time Discrete Gabor Transform for Real signals
 *  ---------------------------------------------------------------
 *
 *
 *
 *
 *
 * @file rtdgtreal.h
 * @author Zdeněk Průša
 * @date 29 Mar 2016
 * @brief Real Time DGTREAL header
 */

#ifndef _RTDGTREAL_H
#define _RTDGTREAL_H

typedef enum
{
    LTFAT_RTDGTPHASE_ZERO,
    LTFAT_RTDGTPHASE_HALFSHIFT
} rtdgt_phasetype;

typedef enum
{
    LTFAT_RTDGT_STREAMOK = 0,
    LTFAT_RTDGT_STREAMOVERFLOW = 1,
    LTFAT_RTDGT_STREAMUNDERFLOW = 2
} rtdgt_streamstate;

typedef enum
{
    LTFAT_FORWARD,
    LTFAT_INVERSE
} ltfat_transformdirection;

#endif /* _RTDGTREAL_H */

typedef struct LTFAT_NAME(rtdgtreal_plan) LTFAT_NAME(rtdgtreal_plan);
// For now, the inverse plan is the same
typedef LTFAT_NAME(rtdgtreal_plan) LTFAT_NAME(rtidgtreal_plan);

/** Create RTDGTREAL plan
 *
 * The function returns NULL if the FFTW plan cannot be crated or there is not enough
 * memory to allocate internal buffers.
 *
 * \param[in]  g      Window
 * \param[in]  gl     Window length
 * \param[in]  M      Number of FFT channels
 * \param[in]  ptype  Phase convention
 *
 * \returns RTDGTREAL plan
 */
LTFAT_EXTERN LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                           const ltfatInt M, const rtdgt_phasetype ptype);

/** Execute RTDGTREAL plan
 * \param[in]  p      RTDGTREAL plan
 * \param[in]  f      Input buffer (gl x W)
 * \param[in]  W      Number of channels
 * \param[out] c      Output DGT coefficients (M2 x W)
 */
LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_execute)(const LTFAT_NAME(rtdgtreal_plan)* p,
                              const LTFAT_REAL* f, const ltfatInt W,
                              LTFAT_COMPLEX* c);

/** Destroy RTDGTREAL plan
 * \param[in]  p      RTDGTREAL plan
 */
LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_done)(LTFAT_NAME(rtdgtreal_plan)* p);

/** Create RTIDGTREAL plan
 *
 * The function returns NULL if the FFTW plan cannot be crated or there is not enough
 * memory to allocate internal buffers.
 *
 * \param[in]  g      Window
 * \param[in]  gl     Window length
 * \param[in]  M      Number of FFT channels
 * \param[in]  ptype  Phase convention
 *
 * \returns RTIDGTREAL plan
 */
LTFAT_EXTERN LTFAT_NAME(rtidgtreal_plan)*
LTFAT_NAME(rtidgtreal_init)(const LTFAT_REAL* g, const ltfatInt gl,
                            const ltfatInt M, const rtdgt_phasetype ptype);

/** Execute RTIDGTREAL plan
 * \param[in]  p      RTDGTREAL plan
 * \param[int] c      Input DGT coefficients (M2 x W)
 * \param[in]  W      Number of channels
 * \param[out] f      Output buffer (gl x W)
 */
LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_execute)(const LTFAT_NAME(rtidgtreal_plan)* p,
                               const LTFAT_COMPLEX* c, const ltfatInt W,
                               LTFAT_REAL* f);

/** Destroy RTIDGTREAL plan
 * \param[in]  p      RTIDGTREAL plan
 */
LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_done)(LTFAT_NAME(rtidgtreal_plan)* p);


LTFAT_NAME(rtdgtreal_plan)*
LTFAT_NAME(rtdgtreal_commoninit)(const LTFAT_REAL* g, const ltfatInt gl,
                                 const ltfatInt M, const rtdgt_phasetype ptype,
                                 const  ltfat_transformdirection tradir);


typedef struct LTFAT_NAME(rtdgtreal_fifo) LTFAT_NAME(rtdgtreal_fifo);

/** Create ring buffer for DGT analysis
 *
 * The ring buffer works as usual when written to, but only constant
 * size (gl) chunks can be read from it and the read pointer is
 * only advanced by a after read.
 *
 * The buffer read and write pointers are initialized such that they
 * reflect the processing delay.
 *
 * \param[in]  fifoLen  Ring buffer size. This should be at least gl + max. expected
 *                      buffer length.
 *                      One more slot is actually allocated for the
 *                      "one slot open" implementation.
 * \param[in]  gl       Window length
 * \param[in]  a        Hop factor
 * \param[in]  Wmax     Maximum number of channels
 *
 * \returns RTDGTREAL_FIFO struct pointer
 */
LTFAT_EXTERN LTFAT_NAME(rtdgtreal_fifo)*
LTFAT_NAME(rtdgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl, ltfatInt a,
                                ltfatInt Wmax);

/** Write bufLen samples to DGT analysis ring buffer
 *
 * The function returns 0 if bufLen is <=0 or the buffer is full.
 * If there is not enough space for all bufLen samples, only available space is used
 * and the number of actually written samples is returned.
 *
 * \param[in]  p        Analysis ring buffer struct
 * \param[in]  buf      Channels to be written.
 * \param[in]  bufLen   Number of samples to be written
 * \param[in]  W        Number of channels
 *
 * \returns Number of samples written
 */
LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_write)(LTFAT_NAME(rtdgtreal_fifo)* p, const LTFAT_REAL** buf,
                                 const ltfatInt bufLen, const ltfatInt W);

/** Read p->gl samples from DGT analysis ring buffer
 *
 * The function attempts to read p->gl samples from the buffer.
 *
 * The function returns 0 if there is less than p->gl samples available.
 *
 * \param[in]   p        Analysis ring buffer struct
 * \param[out]  buf      Output array, it is expected to be able to hold p->gl*p->Wmax samples.
 *
 * \returns Number of samples read
 */
LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_read)(LTFAT_NAME(rtdgtreal_fifo)* p, LTFAT_REAL* buf);

/** Destroy DGT analysis ring buffer
 * \param[in]  p      DGT analysis ring buffer
 */
LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_fifo_done)(LTFAT_NAME(rtdgtreal_fifo)* p);


typedef struct LTFAT_NAME(rtidgtreal_fifo) LTFAT_NAME(rtidgtreal_fifo);

/** Create ring buffer for DGT synthesis
 *
 * The ring buffer behaves as usual when read from, except it sets the read
 * samples to zero.
 * Only chunks of size gl can be written to it and the write pointer is advanced
 * by a. The samples are added to the existing values instead of the usual
 * overwrite.
 *
 * The buffer read and write pointers are both initialized to the same value.
 *
 * \param[in]  fifoLen  Ring buffer size. This should be at least gl + max. expected
 *                      buffer length. (gl+1) more slots are actually allocated
 *                      to accomodate the overlaps.
 * \param[in]  gl       Window length
 * \param[in]  a        Hop factor
 * \param[in]  Wmax     Maximum number of channels
 *
 * \returns RTIDGTREAL_FIFO struct pointer
 */
LTFAT_EXTERN LTFAT_NAME(rtidgtreal_fifo)*
LTFAT_NAME(rtidgtreal_fifo_init)(ltfatInt fifoLen, ltfatInt gl, ltfatInt a, ltfatInt Wmax);

/** Write p->gl samples to DGT synthesis ring buffer
 *
 * The function returns 0 if there is not enough space to write all
 * p->gl samples.
 *
 * \param[in]  p        Synthesis ring buffer struct
 * \param[in]  buf      Samples to be written
 *
 * \returns Number of samples written
 */
LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_write)(LTFAT_NAME(rtidgtreal_fifo)* p, const LTFAT_REAL* buf);

/** Read bufLen samples from DGT analysis ring buffer
 *
 * The function attempts to read bufLen samples from the buffer.
 *
 * \param[in]   p        Analysis ring buffer struct
 * \param[in]   bufLen   Number of samples to be read
 * \param[in]   W        Number of channels
 * \param[out]  buf      Output channels, each channel is expected to be able to
 *                       hold bufLen samples.
 *
 * \returns Number of samples read
 */
LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_read)(LTFAT_NAME(rtidgtreal_fifo)* p,
                                 const ltfatInt bufLen, const ltfatInt W,
                                 LTFAT_REAL** buf);

/** Destroy DGT synthesis ring buffer
 * \param[in]  p      DGT synthesis ring buffer
 */
LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_fifo_done)(LTFAT_NAME(rtidgtreal_fifo)* p);

/** Processor callback signature
 * It is safe to assume that out and in are not aliased.
 *
 * \param[in]  userdata   User defined data
 * \param[in]        in   Input coefficients
 * \param[in]        M2   Length of the arrays; number of unique FFT channels; equals to M/2 + 1
 * \param[out]      out   Output coefficients
 */
typedef void (*LTFAT_NAME(rtdgtreal_processor_callback))(void* userdata,
        const LTFAT_COMPLEX* in, const int M2, const int W, LTFAT_COMPLEX* out);

typedef struct LTFAT_NAME(rtdgtreal_processor) LTFAT_NAME(rtdgtreal_processor);

/** Create DGTREAL processor
 *
 * The processor wraps DGTREAL analysis-modify-synthesis loop suitable for
 * stream of data. The -modify- part is user definable via callback.
 * If the callback is NULL, no coefficient modification occurs.
 *
 * \param[in]        g   Analysis window
 * \param[in]       gd   Synthesis window
 * \param[in]       gl   Length of the windows
 * \param[in]        a   Hop size
 * \param[in]        M   Number of FFT channels
 * \param[in]     Wmax   Maximum number of channels
 * \param[in] callback   Custom function to process the coefficients
 * \param[in] userdata   Custom callback data. Will be passed to the callback.
 *                       Useful for storing state between callback calls.
 *
 * \returns DGTREAL processor struct
 */
LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_init)(const LTFAT_REAL* g, const LTFAT_REAL* gd,
                                     const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                     const ltfatInt Wmax,
                                     LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                     void* userdata);

/** Create DGTREAL processor 
 *
 * This function provides an alternative way of creating DGTREAL processor
 * struct. The function accepts only the analysis window and the synthesis 
 * window is computed internally. 
 * 
 * \param[in]      win   Analysis window
 * \param[in]       gl   Length of the windows
 * \param[in]        a   Hop size
 * \param[in]        M   Number of FFT channels
 * \param[in]     Wmax   Maximum number of channels
 * \param[in] callback   Custom function to process the coefficients
 * \param[in] userdata   Custom callback data. Will be passed to the callback.
 *                       Useful for storing state between callback calls.
 *
 * \returns DGTREAL processor struct
 */
LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_wininit)(LTFAT_FIRWIN win,
                                        const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                        const ltfatInt Wmax,
                                        LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                        void* userdata);

/** Execute DGTREAL processor
 *
 * Output is lagging behind the input by (gl-1) samples.
 * The function can run inplace i.e. in==out.
 *
 * \param[in]      p  DGTREAL processor
 * \param[in]     in  Input channels
 * \param[in]    len  Length of the channels
 * \param[in] chanNo  Number of channels
 * \param[out]   out  Output frame
 *
 * \returns Error status
 */
LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_processor_execute)(LTFAT_NAME(rtdgtreal_processor)* p,
                                        const LTFAT_REAL** in,
                                        const ltfatInt len, const ltfatInt chanNo,
                                        LTFAT_REAL** out);

/** Set DGTREAL processor callback 
 *
 * Function replaces the callback in the struct. This is not thread safe.
 * Only call this if there is no chance that the execute function is called
 * simultaneously in a different thread.
 *
 * \param[in/out]        p   DGTREAL processor
 * \param[in]     callback   Custom function to process the coefficients
 * \param[in]     userdata   Custom callback data. Will be passed to the callback.
 *                           Useful for storing state between callback calls.
 *
 * \returns DGTREAL processor struct
 */
LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_processor_setcallback)(LTFAT_NAME(rtdgtreal_processor)* p,
                             LTFAT_NAME(rtdgtreal_processor_callback) callback,
                             void* userdata);

/** Destroy DGTREAL processor
 * \param[in]  p      DGTREAL processor
 */
LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_processor_done)(LTFAT_NAME(rtdgtreal_processor)* p);

/** Default processor callback
 * \param[in]  userdata   User defined data (unused)
 * \param[in]        in   Input coefficients
 * \param[in]        M2   Length of the arrays; number of unique FFT channels; equals to M/2 + 1
 * \param[int]        W   Number of channels
 * \param[out]      out   Output coefficients
 */
LTFAT_EXTERN void
LTFAT_NAME(default_rtdgtreal_processor_callback)(void* userdata, const LTFAT_COMPLEX* in,
        const int M2, const int W, LTFAT_COMPLEX* out);

/** @}*/

