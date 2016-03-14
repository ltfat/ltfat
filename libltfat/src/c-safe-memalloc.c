#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN
void* ltfat_malloc (size_t n)
{
    void *outp;
    outp = LTFAT_FFTW(malloc)(n);
    if (outp == NULL)
    {
        puts("ltfat_malloc failed.");
        exit(1);
    }

    return outp;
}

LTFAT_EXTERN
void* ltfat_realloc (void *ptr, size_t n)
{
    void *outp;
    // DOES NOT PRODUCE MEMORY ALIGNED POINTER
    // outp = realloc(ptr, n);
    outp = LTFAT_FFTW(malloc)(n);

    if (outp == NULL)
    {
        puts("ltfat_realloc failed.");
        exit(1);
    }

    if(ptr!=NULL)
    {
        ltfat_free(ptr);
    }

    return outp;
}

void* ltfat_realloc_and_copy (void *ptr, size_t nold, size_t nnew)
{
    if (ptr == NULL)
    {
        puts("Null pointer.");
        exit(1);
    }

    void *outp;

    outp = LTFAT_FFTW(malloc)(nnew);

    if (outp == NULL)
    {
        puts("ltfat_realloc_and_copy failed.");
        exit(1);
    }

    memcpy(outp,ptr,nold<nnew?nold:nnew);

    ltfat_free(ptr);

    return outp;
}

LTFAT_EXTERN
void* ltfat_calloc (size_t nmemb, size_t size)
{
    void *outp;
    // DOES NOT PRODUCE MEMORY ALIGNED POINTER
    // outp = calloc(nmemb, size);

    // workaround
    outp = LTFAT_FFTW(malloc)(nmemb*size);

    if (outp == NULL)
    {
        puts("ltfat_calloc failed.");
        exit(1);
    }
    // workaround
    memset(outp,0,nmemb*size);

    return outp;
}

LTFAT_EXTERN
void ltfat_free(const void *ptr)
{
    LTFAT_FFTW(free)((void*)ptr);
}

void ltfat_safefree(const void *ptr)
{
    if(ptr!=NULL)
        ltfat_free((void *)ptr);
}
