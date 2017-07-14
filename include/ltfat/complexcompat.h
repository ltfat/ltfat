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
#include <fftw3.h>
#endif

#if defined(__cplusplus) && !defined(LTFAT_CINTERFACE)
    typedef std::complex<double> ltfat_complex_d;
    typedef std::complex<float> ltfat_complex_s;
#else
#   if defined(_Complex_I) && defined(complex) && defined(I)
        typedef complex double ltfat_complex_d;
        typedef complex float ltfat_complex_s;
#   else
        typedef double ltfat_complex_d[2];
        typedef float ltfat_complex_s[2] ;
#   endif
#endif

// #if defined(__cplusplus)
// #   define ltfat_real(x) std::real(x)
// #   define ltfat_imag(x) std::imag(x)
// #   define ltfat_abs(x) std::abs(x)
// #   define ltfat_arg(x) std::arg(x)
// #else
// #   define ltfat_complex_d(r,i) ((float)(r) + ((float)(i))*I)
// #   define ltfat_complex_s(r,i) ((double)(r) + ((double)(i))*I)
// #   define ltfat_real(x) creal(x)
// #   define ltfat_imag(x) cimag(x)
// #   define ltfat_abs(x) fabs(x)
// #   define ltfat_arg(x) carg(x)
// #endif

#endif
