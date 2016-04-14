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

typedef struct
{
    const LTFAT_REAL* g; //!< Window
    const ltfatInt gl; //!< Window length
    const ltfatInt M; //!< Number of FFT channels
    const rtdgt_phasetype ptype; //!< Phase convention
    LTFAT_REAL* fftBuf; //!< Internal buffer
    const ltfatInt fftBufLen; //!< Internal buffer length
    LTFAT_FFTW(plan) pfft; //!< FFTW plan
} LTFAT_NAME(rtdgtreal_plan);

// For now, the inverse plan is teh same
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
 * \param[out] c      Output DGT coefficients
 * \param[in]  ptype  Phase convention
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
 * \param[in]  f      Input buffer (gl x W)
 * \param[in]  W      Number of channels
 * \param[out] c      Output DGT coefficients
 * \param[in]  ptype  Phase convention
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


typedef struct
{
    ltfatInt gl;
    ltfatInt a;
    LTFAT_REAL* buf;
    ltfatInt bufLen;
    ltfatInt readIdx;
    ltfatInt writeIdx;
} LTFAT_NAME(rtdgtreal_fifo);


LTFAT_EXTERN LTFAT_NAME(rtdgtreal_fifo)*
LTFAT_NAME(rtdgtreal_fifo_init)(ltfatInt bufLen, ltfatInt gl, ltfatInt a);


LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_write)(LTFAT_NAME(rtdgtreal_fifo)* p, const ltfatInt bufLen, const LTFAT_REAL* buf);

LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_fifo_read)(LTFAT_NAME(rtdgtreal_fifo)* p, LTFAT_REAL* buf);


LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_fifo_done)(LTFAT_NAME(rtdgtreal_fifo)* p);


typedef struct
{
    ltfatInt gl;
    ltfatInt a;
    LTFAT_REAL* buf;
    ltfatInt bufLen;
    ltfatInt readIdx;
    ltfatInt writeIdx;
} LTFAT_NAME(rtidgtreal_fifo);

LTFAT_EXTERN LTFAT_NAME(rtidgtreal_fifo)*
LTFAT_NAME(rtidgtreal_fifo_init)(ltfatInt bufLen, ltfatInt gl, ltfatInt a);

LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_write)(LTFAT_NAME(rtidgtreal_fifo)* p, const LTFAT_REAL* buf);

LTFAT_EXTERN int
LTFAT_NAME(rtidgtreal_fifo_read)(LTFAT_NAME(rtidgtreal_fifo)* p, const ltfatInt bufLen, LTFAT_REAL* buf);


LTFAT_EXTERN void
LTFAT_NAME(rtidgtreal_fifo_done)(LTFAT_NAME(rtidgtreal_fifo)* p);


typedef void (*LTFAT_NAME(rtdgtreal_processor_callback))(void* userdata, const LTFAT_COMPLEX* in, const ltfatInt M2, LTFAT_COMPLEX* out);

typedef struct
{
    LTFAT_NAME(rtdgtreal_processor_callback) processorCallback;
    void* userdata;
    LTFAT_NAME(rtdgtreal_fifo)* fwdfifo;
    LTFAT_NAME(rtidgtreal_fifo)* backfifo;
    LTFAT_NAME(rtdgtreal_plan)* fwdplan;
    LTFAT_NAME(rtidgtreal_plan)* backplan;
    LTFAT_REAL* buf;
    LTFAT_COMPLEX* fftbufIn;
    LTFAT_COMPLEX* fftbufOut;
} LTFAT_NAME(rtdgtreal_processor);

LTFAT_EXTERN LTFAT_NAME(rtdgtreal_processor)*
LTFAT_NAME(rtdgtreal_processor_init)(const LTFAT_REAL* g, const LTFAT_REAL* gd,
                                     const ltfatInt gl, const ltfatInt a, const ltfatInt M,
                                     LTFAT_NAME(rtdgtreal_processor_callback) callback,
                                     void* userdata);

// LTFAT_EXTERN int
// LTFAT_NAME(rtdgtreal_processor_registercallback)(LTFAT_NAME(rtdgtreal_processor)* p,
//         LTFAT_NAME(rtdgtreal_processor_callback) callback, void* userdata);

LTFAT_EXTERN int
LTFAT_NAME(rtdgtreal_processor_execute)(LTFAT_NAME(rtdgtreal_processor)* p,
                                        const LTFAT_REAL* in, const ltfatInt len, LTFAT_REAL* out);

LTFAT_EXTERN void
LTFAT_NAME(rtdgtreal_processor_done)(LTFAT_NAME(rtdgtreal_processor)* p);

LTFAT_EXTERN void
LTFAT_NAME(default_rtdgtreal_processor_callback)(void* userdata, const LTFAT_COMPLEX* in, const ltfatInt M2, LTFAT_COMPLEX* out);

/** @}*/
