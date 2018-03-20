#ifndef _LTFAT_COMPLEXCOMPAT_H
#define _LTFAT_COMPLEXCOMPAT_H

#ifndef NOSYSTEMHEADERS

#ifdef __cplusplus
// C++ complex header
// fftw3.h will define:
// typedef double fftw_complex[2]
#include <complex>
#include <cmath>
using namespace std;
#else
// C99 complex header
// fftw3.h will define:
// typedef double _Complex fftw_complex
// #include <complex.h>
// #include <math.h>
#include <tgmath.h>
#endif

// Must be included after complex.h
// #include "ltfat/thirdparty/fftw3.h"
#endif

#if defined(__cplusplus) && !defined(LTFAT_CINTERFACE)
    typedef std::complex<double> ltfat_complex_d;
    typedef std::complex<float> ltfat_complex_s;
#else
#if !defined(LTFAT_CINTERFACE) && defined(_Complex_I) && defined(complex) && defined(I)
    typedef double _Complex ltfat_complex_d;
    typedef float  _Complex ltfat_complex_s;
#else
    typedef double ltfat_complex_d[2];
    typedef float  ltfat_complex_s[2];
#endif
#endif

#endif
