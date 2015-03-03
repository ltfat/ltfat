/**
*
*/

#ifdef LTFAT_COMPLEX
#undef LTFAT_COMPLEX
#endif
#ifdef LTFAT_REAL
#undef LTFAT_REAL
#endif
#ifdef LTFAT_TYPE
#undef LTFAT_TYPE
#endif
#ifdef LTFAT_NAME
#undef LTFAT_NAME
#endif
#ifdef LTFAT_NAME_REAL
#undef LTFAT_NAME_REAL
#endif
#ifdef LTFAT_NAME_COMPLEX
#undef LTFAT_NAME_COMPLEX
#endif
#ifdef LTFAT_FFTW
#undef LTFAT_FFTW
#endif

#ifdef LTFAT_MX_CLASSID
#undef LTFAT_MX_CLASSID
#endif

#ifdef LTFAT_MX_COMPLEXITY
#undef LTFAT_MX_COMPLEXITY
#endif

#ifdef LTFAT_COMPLEXH
#undef LTFAT_COMPLEXH
#endif

#ifdef LTFAT_COMPLEXH
#undef LTFAT_COMPLEXH
#endif

#ifdef LTFAT_DOUBLE
#  define LTFAT_REAL double
#  define LTFAT_COMPLEX fftw_complex
#  define LTFAT_FFTW(name) fftw_ ## name
#  define LTFAT_NAME_REAL(name) LTFAT_NAME_DOUBLE(name)
#  define LTFAT_NAME_COMPLEX(name) LTFAT_NAME_COMPLEXDOUBLE(name)
#  define LTFAT_COMPLEXH(name) name
#  define LTFAT_MX_CLASSID mxDOUBLE_CLASS
#  if defined(LTFAT_COMPLEXTYPE)
#    define LTFAT_TYPE LTFAT_COMPLEX
#    define LTFAT_NAME(name) LTFAT_NAME_COMPLEXDOUBLE(name)
#    define LTFAT_MX_COMPLEXITY mxCOMPLEX
#  else
#    define LTFAT_TYPE LTFAT_REAL
#    define LTFAT_NAME(name) LTFAT_NAME_DOUBLE(name)
#    define LTFAT_MX_COMPLEXITY mxREAL
#  endif
#endif

#ifdef LTFAT_SINGLE
#define LTFAT_REAL float
#define LTFAT_COMPLEX fftwf_complex
#define LTFAT_MX_CLASSID mxSINGLE_CLASS
#define LTFAT_NAME_REAL(name) LTFAT_NAME_SINGLE(name)
#define LTFAT_NAME_COMPLEX(name) LTFAT_NAME_COMPLEXSINGLE(name)
#define LTFAT_FFTW(name) fftwf_ ## name
#define LTFAT_COMPLEXH(name) name ## f
#  if defined(LTFAT_COMPLEXTYPE)
#    define LTFAT_TYPE LTFAT_COMPLEX
#    define LTFAT_NAME(name) LTFAT_NAME_COMPLEXSINGLE(name)
#    define LTFAT_MX_COMPLEXITY mxCOMPLEX
#  else
#    define LTFAT_TYPE LTFAT_REAL
#    define LTFAT_NAME(name) LTFAT_NAME_SINGLE(name)
#    define LTFAT_MX_COMPLEXITY mxREAL
#  endif
#endif


