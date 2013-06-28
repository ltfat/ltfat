AC_PREREQ([2.62])
AC_INIT([LTFAT], [1.4.0], [http://ltfat.sourceforge.net/contact.php], [ltfat])

AC_PROG_CC
AC_PROG_CXX

AC_CONFIG_MACRO_DIR([m4/])

dnl Checks for standard programs.
AC_PROG_F77
AC_PROG_LIBTOOL

dnl Check for BLAS libraries
sinclude(ax_blas.m4)
AX_BLAS
if test "$ax_blas_ok" = "no"; then
   AC_MSG_ERROR([Cannot find BLAS libraries])
fi

dnl Check for LAPACK libraries
sinclude(ax_lapack.m4)
AX_LAPACK
if test "$ax_lapack_ok" = "no"; then
   AC_MSG_ERROR([Cannot find LAPACK libraries])
fi

# Put in check for MKOCTFILE

AC_CHECK_HEADERS([fftw3.h],[],[AC_MSG_ERROR([fftw was not found])])
# Put in check for blas
# Put in check for lapack
# Put in check for portaudio

AC_OUTPUT
