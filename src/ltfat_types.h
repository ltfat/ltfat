  //#include "fftw3.h"

  #ifdef LTFAT_COMPLEX
     #undef LTFAT_COMPLEX
  #endif
  #ifdef LTFAT_REAL
     #undef LTFAT_REAL
  #endif
  #ifdef LTFAT_NAME
     #undef LTFAT_NAME
  #endif
  #ifdef LTFAT_FFTW
     #undef LTFAT_FFTW
  #endif

  #ifdef LTFAT_MX_CLASSID
     #undef LTFAT_MX_CLASSID
  #endif


#ifdef LTFAT_DOUBLE
#define LTFAT_COMPLEX fftw_complex
#define LTFAT_REAL double
#define LTFAT_NAME(name) name
#define LTFAT_FFTW(name) fftw_ ## name
#define LTFAT_MX_CLASSID mxDOUBLE_CLASS
#endif

#ifdef LTFAT_SINGLE
#define LTFAT_COMPLEX fftwf_complex
#define LTFAT_REAL float
#define LTFAT_NAME(name) s ## name
#define LTFAT_FFTW(name) fftwf_ ## name
#define LTFAT_MX_CLASSID mxSINGLE_CLASS
#endif
