#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN
void* ltfat_malloc (size_t n)
{
    void* outp;
    outp = LTFAT_FFTW(malloc)(n);
    return outp;
}

LTFAT_EXTERN
void* ltfat_realloc (void* ptr, size_t n)
{
    void* outp;
    // DOES NOT PRODUCE MEMORY ALIGNED POINTER
    // outp = realloc(ptr, n);
    outp = LTFAT_FFTW(malloc)(n);

    if (!outp)
        return NULL;

    if (ptr)
        ltfat_free(ptr);

    return outp;
}

void* ltfat_realloc_and_copy (void* ptr, size_t nold, size_t nnew)
{
    void* outp = LTFAT_FFTW(malloc)(nnew);

    if (!outp)
        return NULL;

    if (ptr)
    {
        memcpy(outp, ptr, nold < nnew ? nold : nnew);
        ltfat_free(ptr);
    }

    return outp;
}

LTFAT_EXTERN
void* ltfat_calloc (size_t nmemb, size_t size)
{
    void* outp = LTFAT_FFTW(malloc)(nmemb * size);

    if (!outp)
        return NULL;

    memset(outp, 0, nmemb * size);

    return outp;
}

LTFAT_EXTERN
void ltfat_free(const void* ptr)
{
    LTFAT_FFTW(free)((void*)ptr);
}

void ltfat_safefree(const void* ptr)
{
    if (ptr)
        ltfat_free((void*)ptr);
}
