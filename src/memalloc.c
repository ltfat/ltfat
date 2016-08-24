#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/memalloc.h"
/* #include "ltfat/macros.h" */

void* (*ltfat_custom_malloc)(size_t) = NULL;
void (*ltfat_custom_free)(void*) = NULL;

ltfat_memory_handler_t
ltfat_set_memory_handler (ltfat_memory_handler_t new_handler)
{
    ltfat_memory_handler_t retVal = { ltfat_custom_malloc, ltfat_custom_free };

    ltfat_custom_malloc = new_handler.malloc;
    ltfat_custom_free = new_handler.free;
    return retVal;
}

LTFAT_API void*
ltfat_malloc (size_t n)
{
    void* outp;

    if (ltfat_custom_malloc)
        outp = (*ltfat_custom_malloc)(n);
    else
        outp = LTFAT_FFTW(malloc)(n);
    return outp;
}

LTFAT_API void*
ltfat_realloc (void* ptr, size_t nold, size_t nnew)
{
    void* outp = ltfat_malloc(nnew);

    if (!outp)
        return NULL;

    if (ptr)
    {
        memcpy(outp, ptr, nold < nnew ? nold : nnew);
        ltfat_free(ptr);
    }

    return outp;
}

LTFAT_API void*
ltfat_calloc (size_t nmemb, size_t size)
{
    void* outp = ltfat_malloc(nmemb * size);

    if (!outp)
        return NULL;

    memset(outp, 0, nmemb * size);

    return outp;
}

LTFAT_API void
ltfat_free(const void* ptr)
{
    if (ltfat_custom_free)
        (*ltfat_custom_free)((void*)ptr);
    else
        LTFAT_FFTW(free)((void*)ptr);
}

LTFAT_API void
ltfat_safefree(const void* ptr)
{
    if (ptr)
        ltfat_free((void*)ptr);
}
