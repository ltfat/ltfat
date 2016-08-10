#ifndef _BASICMACROS_H
#define _BASICMACROS_H

/* Handle Windows DLL files */
/* defined by Makefile when compiling LTFAT */
#if defined(DLL_EXPORT_SYMBOLS) && ((defined(_WIN32) || defined(__WIN32__)))
#  define LTFAT_EXTERN extern __declspec(dllexport)
#else
#  define LTFAT_EXTERN extern
#endif

#define LTFAT_MAKENAME(name,suffix) ltfat ## _ ## name ## suffix
#define LTFAT_NAME_DOUBLE(name) LTFAT_MAKENAME(name,_d)
#define LTFAT_NAME_SINGLE(name) LTFAT_MAKENAME(name,_s)
#define LTFAT_NAME_COMPLEXDOUBLE(name) LTFAT_MAKENAME(name,_dc)
#define LTFAT_NAME_COMPLEXSINGLE(name) LTFAT_MAKENAME(name,_sc)


#endif
