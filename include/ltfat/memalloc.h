/** \defgroup memalloc Memory allocation
 *
 * Internally, heap memory allocation in libltfat is done using functions 
 * in this module.
 *
 * By default,
 * <a href="http://www.fftw.org/doc/Memory-Allocation.html">malloc/free</a>
 * functions from FFTW3 are used. FFTW malloc allocates memory properly aligned.
 *
 * A custom malloc/free can be registered using ltfat_set_memory_handler.
 *
 * \note According to the FFTW documentation, the allocated memory __is not required__
 * to be aligned, but it is __strongly recommended__ to do so.
 *
 * \addtogroup memalloc
 * @{
 */
#ifndef _LTFAT_MEMALLOC_H
#define _LTFAT_MEMALLOC_H


#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
    void* (*malloc) (size_t n);
    void (*free) (void*);
} ltfat_memory_handler_t;

/** Set custom malloc/free functions
 \returns Old malloc/free
 */
ltfat_memory_handler_t
ltfat_set_memory_handler (ltfat_memory_handler_t new_handler);

LTFAT_EXTERN
void* ltfat_malloc (size_t n);

LTFAT_EXTERN
void* ltfat_calloc (size_t nmemb, size_t size);

LTFAT_EXTERN
void* ltfat_realloc (void *ptr, size_t nold, size_t nnew);

LTFAT_EXTERN
void  ltfat_free(const void *ptr);

/** @} */

LTFAT_EXTERN
void  ltfat_safefree(const void *ptr);

#ifdef __cplusplus
}
#endif



#endif
