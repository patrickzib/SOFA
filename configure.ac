AC_INIT([ParIS Parallel index of sequence], [1.3], [botao.peng@parisdescartes.fr],
             [libads])
AC_PREREQ([2.59])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])
AC_CONFIG_HEADERS([config.h])
AM_PROG_AR
AC_PROG_RANLIB
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

AC_CHECK_HEADER(values.h)
AC_CHECK_HEADER(fftw3.h)
AC_CANONICAL_HOST
CANONICAL_HOST=$host
case $host in
	*-*-linux*)
		AC_CHECK_HEADER(numa.h)
		LIB64='64'
		;;
	*-*-darwin*)
		echo "darwin detected, skipping numa.h dep, setting brew environment"
		BREW_ROOT=`brew config|grep HOMEBREW_PREFIX|awk '{print $2}'`
		APPLE_OMP='-I${BREW_ROOT}/opt/libomp/include'
		APPLE_FFTW='-I${BREW_ROOT}/opt/fftw/include'
		APPLE_NOWARN_CONV='-Wno-int-conversion'
		LDFLAGS='${LD_FLAGS} -L${BREW_ROOT}/opt/fftw/lib'
		case $host_cpu in
			arm)
				BUILD_ARCH='armv8.2-a'
				;;
		esac
		;;
esac
case $host_cpu in
	x86_64)
		echo "x86_64-cpu"
		BUILD_ARCH='haswell'
		NUMA='-lnuma'
		;;
	aarch64)
		echo "arm64-cpu"
		BUILD_ARCH='armv8.2-a'
		;;
esac
AC_SUBST([LIB64])
AC_SUBST([BUILD_ARCH])
AC_SUBST([NUMA])
AC_SUBST([BREW_ROOT])
AC_SUBST([APPLE_OMP])
AC_SUBST([APPLE_FFTW])
AC_SUBST([APPLE_NOWARN_CONV])
AC_SUBST([LDFLAGS])

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


	AC_CONFIG_FILES([Makefile])
	AC_OUTPUT
