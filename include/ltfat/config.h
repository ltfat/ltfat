#ifndef CONFIG_H
#define CONFIG_H 1

#define LTFAT_MAKENAME(name,type,comp) name ## _ ## comp ## type
#define LTFAT_NAME_DOUBLE(name) LTFAT_MAKENAME(name,d,)
#define LTFAT_NAME_SINGLE(name) LTFAT_MAKENAME(name,s,)
#define LTFAT_NAME_COMPLEXDOUBLE(name) LTFAT_MAKENAME(name,d,c)
#define LTFAT_NAME_COMPLEXSINGLE(name) LTFAT_MAKENAME(name,s,c)

/* Handle Windows DLL files */
/* defined by Makefile when compiling LTFAT */
#if defined(DLL_EXPORT_SYMBOLS) && ((defined(_WIN32) || defined(__WIN32__)))
#  define LTFAT_EXTERN extern __declspec(dllexport)
#else
#  define LTFAT_EXTERN extern
#endif

/* On WinXP, gcc defines __WIN32__ */
/* On Linux, gcc defines __linux__ */

#endif /* CONFIG_H */
