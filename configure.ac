AC_INIT([ParIS Parallel index of sequence], [1.3], [botao.peng@parisdescartes.fr],
             [libads])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])
AC_CONFIG_HEADERS([config.h])
AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR
AM_PROG_CC_C_O
AC_CHECK_HEADER(values.h)

if test "$ac_cv_header_values_h" == yes
then
	AC_DEFINE(VALUES, [], [include values.h])
fi

AC_ARG_WITH(no_benchmarking, [  --no-benchmarking            disable stuff])
if test "$no_benchmarking" == "yes"
then
	AC_MSG_NOTICE([You have disabled benchmarking!])
else
	AC_DEFINE(BENCHMARK, [], [Create a benchmark version])
fi

AC_ARG_WITH(with_debug, [  --with-debug            enable stuff])
if test "$with_debug" == "yes"
then
	AC_MSG_NOTICE([You have enabled debugging info!])
	AC_DEFINE(DEBUG, [], [Enable debugging info])
fi

AC_ARG_WITH(with_clustering, [  --with-clustering            enable stuff])
if test "$with_clustering" == "yes"
then
	AC_MSG_NOTICE([You have enabled input file clustering!])
	AC_DEFINE(CLUSTERED, [], [Enable clustering])
fi

AC_ARG_WITH(verbose, [  --verbose            enable stuff])
if test "$verbose" == "yes"
then
	AC_MSG_NOTICE([You have enabled verbose output!])
	AC_DEFINE(VERBOSE, [2], [Enable verbose output])
else
	AC_DEFINE(VERBOSE, [0], [Disable verbose output])
fi

AC_ARG_ENABLE([simd],
    [AS_HELP_STRING([--enable-simd@<:@=auto@:>@],
        [enable AVX2 SIMD optimizations (default=auto; use 'no' to disable)])],
    [],
    [enable_simd=auto])

SIMD_CFLAGS=""
ads_have_avx2=0

if test "x$enable_simd" != "xno"; then
    AC_MSG_CHECKING([whether the compiler supports AVX2 intrinsics])
    simd_save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS -mavx2"
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[__m256 a = _mm256_set1_ps(1.0f);
                            __m256 b = _mm256_add_ps(a, a);
                            return (int)_mm_cvtss_f32(_mm256_castps256_ps128(b));]])],
        [ads_have_avx2=1],
        [ads_have_avx2=0])
    CFLAGS="$simd_save_CFLAGS"
    if test $ads_have_avx2 -eq 1; then
        AC_MSG_RESULT([yes])
        SIMD_CFLAGS="-mavx -mavx2 -msse3"
        AC_DEFINE([ADS_HAVE_AVX2], [1], [Define if AVX2 SIMD is available])
    else
        AC_MSG_RESULT([no])
        if test "x$enable_simd" = "xyes"; then
            AC_MSG_ERROR([--enable-simd was requested, but AVX2 intrinsics are not supported by this compiler])
        fi
    fi
fi

if test $ads_have_avx2 -ne 1; then
    AC_DEFINE([ADS_HAVE_AVX2], [0], [Define if AVX2 SIMD is available])
fi

AC_SUBST([SIMD_CFLAGS])

PKG_PROG_PKG_CONFIG
AC_ARG_WITH([fftw],
    [AS_HELP_STRING([--with-fftw=PREFIX], [Prefix where FFTW3 (float) is installed])],
    [FFTW_CPPFLAGS="-I$with_fftw/include"
     FFTW_LDFLAGS="-L$with_fftw/lib"],
    [FFTW_CPPFLAGS=""
     FFTW_LDFLAGS=""])

save_CPPFLAGS="$CPPFLAGS"
save_LDFLAGS="$LDFLAGS"
save_LIBS="$LIBS"

CPPFLAGS="$CPPFLAGS $FFTW_CPPFLAGS"
LDFLAGS="$LDFLAGS $FFTW_LDFLAGS"

fftw_found=no

if test "x$PKG_CONFIG" != "x" && $PKG_CONFIG --exists fftw3f; then
    FFTW_CFLAGS="$($PKG_CONFIG --cflags fftw3f)"
    FFTW_LIBS="$($PKG_CONFIG --libs fftw3f)"
    fftw_found=yes
else
    AC_CHECK_HEADER([fftw3.h], [],
        [AC_MSG_ERROR([FFTW3 header fftw3.h not found; install FFTW 3.x or use --with-fftw=PREFIX])])
    AC_CHECK_LIB([fftw3f], [fftwf_plan_dft_1d],
        [FFTW_LIBS="-lfftw3f"; fftw_found=yes],
        [AC_MSG_ERROR([FFTW3 library (fftw3f) not found; install FFTW 3.x or use --with-fftw=PREFIX])])
fi

CPPFLAGS="$save_CPPFLAGS"
LDFLAGS="$save_LDFLAGS"
LIBS="$save_LIBS"

AC_SUBST([FFTW_CFLAGS])
AC_SUBST([FFTW_LIBS])

AC_MSG_CHECKING([for pthread_barrier_t])
barrier_save_LIBS="$LIBS"
LIBS="$LIBS -lpthread"
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <pthread.h>]],
                     [[pthread_barrier_t barrier;
                        pthread_barrier_init(&barrier, NULL, 2);
                        pthread_barrier_wait(&barrier);
                        pthread_barrier_destroy(&barrier);
                        return 0;]])],
    [AC_MSG_RESULT([yes])
     AC_DEFINE([HAVE_PTHREAD_BARRIER], [1], [Define if pthread_barrier_t is available])],
    [AC_MSG_RESULT([no])
     AC_DEFINE([HAVE_PTHREAD_BARRIER], [0], [Define if pthread_barrier_t is available])])
LIBS="$barrier_save_LIBS"

	AC_CONFIG_FILES([Makefile])
	AC_OUTPUT
