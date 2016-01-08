AC_PREREQ([2.62])
AC_INIT([LTFAT], [0.0.1], [http://github.com/ltfat/ltfat/issues], [ltfat])

AC_PROG_CC
AC_PROG_CC_C99

#AC_CONFIG_MACRO_DIR([m4/])
m4_include([m4/ax_blas.m4])
m4_include([m4/ax_lapack.m4])

# The checks for BLAS and Lapack have been disabled because they
# provoke an error about a missing install-sh

#dnl Check for BLAS libraries
#sinclude(ax_blas.m4)
#AX_BLAS
#if test "$ax_blas_ok" = "no"; then
#   AC_MSG_ERROR([Cannot find BLAS libraries])
#fi

#dnl Check for LAPACK libraries
#sinclude(ax_lapack.m4)
#AX_LAPACK
#if test "$ax_lapack_ok" = "no"; then
#   AC_MSG_ERROR([Cannot find LAPACK libraries])
#fi

AC_CHECK_HEADERS([fftw3.h],[],[AC_MSG_ERROR([fftw was not found])])

AC_OUTPUT(Makefile)
