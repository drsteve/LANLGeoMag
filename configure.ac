#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_INIT([LanlGeoMag], [1.5.13], [mghenderson@lanl.gov])

AC_PREREQ([2.59])
#THE_PKG_RELEASE=0
#PKG_RELEASE=$THE_PKG_RELEASE
#AC_SUBST([PKG_RELEASE])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AC_LANG([C])
#AC_DEFINE([LGM_EOP_DATA_DIR],["/usr/local/share/LanlGeoMag/EopData"],["Define directory where Eop data can be found"])
#AH_VERBATIM([PDDH], [#define LGM_EOP_DATA_DIR "@datadir@/LanlGeoMag/EopData"])
#AC_DEFINE([PDD], ["PDD @prefix@/data"], ["dude"])
AC_DEFINE([LGM_EOP_DATA_DIR],[(DATADIR "/LanlGeoMag/EopData")],[Define directory where Eop data can be found])
AC_DEFINE([LGM_INDEX_DATA_DIR],[(DATADIR "/LanlGeoMag/Data")],[Define directory where Geophysical Index Data (e.g. Kp, AE, Dst, etc. data can be found])
AC_DEFINE([LGM_TS07_DATA_DIR],[(DATADIR "/LanlGeoMag/Data")],[Define directory where TS07 data files are located])
AC_CONFIG_SRCDIR([libLanlGeoMag/quicksort.c])
AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL


AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

#AM_CONDITIONAL([HAVE_DOXYGEN],
#[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Doxyfile])])
AM_CPPFLAGS=-DDATADIR='$(datadir)'

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(LanlGeoMag, Doxyfile)



AM_CONDITIONAL(LINUX, [test -n "`echo ${host_os} | grep -i linux`"])
AM_CONDITIONAL(APPLE, [test -n "`echo ${host_os} | grep -i apple`"])

DAILY_CRON_DIR=""
case "$host_os" in
linux*) 
    AM_CONDITIONAL(LINUX, true)
    DAILY_CRON_DIR="/etc/cron.daily"
    #AC_DEFINE(LINUX,1,[This is a Linux system])
    #AC_DEFINE_UNQUOTED(LINUX,1,[This is a LINUX system])
    #AM_CONDITIONAL(HAVE_GNU_LD, true)    
    #AM_CONDITIONAL(SOMEBSD, false)
    #XTRA_CPPLIBS=-lstdc++
    #LIBEXT=.so
    ;;
apple*)
    AM_CONDITIONAL(APPLE, true)
    DAILY_CRON_DIR="/etc/periodic/daily"
    #AC_DEFINE_UNQUOTED(APPLE,1,[This is an Apple system])
    #AM_CONDITIONAL(HAVE_GNU_LD, false)
    #AM_CONDITIONAL(SOMEBSD, false)
    #CFLAGS="-fno-common $CFLAGS"
    #LIBEXT=.so
    ;;
esac




dnl Example of default-enabled feature
AC_ARG_ENABLE([eopcron], AS_HELP_STRING([--enable-eopcron], [Enable running GetEopDataFiles via cron]))

ENABLE_EOP_CRON=0
AS_IF([test "x$enable_eopcron" == "xyes"], [
    ENABLE_EOP_CRON=1
])
AC_SUBST([ENABLE_EOP_CRON])
AC_SUBST([DAILY_CRON_DIR])

AC_ARG_ENABLE([tools], AS_HELP_STRING([--enable-tools], [Install some tools which use LanlGeoMag (default: disable)]))
AS_IF([test "x$enable_tools" == "xyes"],
    [AM_CONDITIONAL([ENABLE_TOOLS], true)],
    [AM_CONDITIONAL([ENABLE_TOOLS], false)]
)

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
dnl Need to directly specify SPICE location for MagEphem
dnl Can change this later to go hunting for them, but standard ld search path won't find it.
AC_ARG_WITH([cspice],
  [AS_HELP_STRING([--with-cspice=DIR], [Path to NAIF/cSPICE (contains include, lib).])],
  [cspice_CFLAGS=-I${withval}/include; cspice_LIBS=-L${withval}/lib/cspice.a] )
AC_ARG_WITH([cspice-includedir],
  [AS_HELP_STRING([--with-cspice-includedir=DIR], [Path to NAIF/cSPICE includes (contains SpiceUsr.h). Required if installing tools.])],
  [cspice_CFLAGS="-I$withval"]
)
AC_ARG_WITH([cspice-libdir],
  [AS_HELP_STRING([--with-cspice-libdir=DIR], [Path to NAIF/cSPICE library (contains cspice.a). Required if installing tools.])],
  [cspice_LIBS="${withval}/cspice.a"]
)
if test "x$enable_tools" == "xyes" -a \( "x$cspice_CFLAGS" == "x" -o "x$cspice_LIBS" == "x" \); then
   AC_MSG_FAILURE([--enable-tools requires either --with-cspice or both --with-cspice-libdir and --with-cspice-includedir])
fi
AC_SUBST(cspice_CFLAGS)
AC_SUBST(cspice_LIBS)

dnl Get flags for embedded perl and make sure it works
dnl Save old CFLAGS/LDFLAGS so can restore after the check
AC_CHECK_PROG([HAVE_PERL], [perl], true, false)
if test x"${HAVE_PERL}" == x"false"; then
    AC_MSG_FAILURE([Unable to run perl])
fi
PERL_CFLAGS=`perl -MExtUtils::Embed -e ccopts`
PERL_LDFLAGS=`perl -MExtUtils::Embed -e ldopts`
OLD_CFLAGS=${CFLAGS}
OLD_LDFLAGS=${LDFLAGS}
OLD_LIBS=${LIBS}
CFLAGS=${CFLAGS}\ ${PERL_CFLAGS}
LDFLAGS=${LDFLAGS}\ ${PERL_LDFLAGS}
AC_CHECK_LIB([perl], [perl_alloc], [AC_DEFINE([HAVE_LIBPERL], [1], [libperl found])], [AC_MSG_ERROR([No libperl found, are perl development files installed?])])
CFLAGS=${OLD_CFLAGS}
LDFLAGS=${OLD_LDFLAGS}
LIBS=${OLD_LIBS}
AC_SUBST(PERL_CFLAGS)
AC_SUBST(PERL_LDFLAGS)

dnl Check for gsl, http://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html
AC_CHECK_LIB([m], [cos], ,
    [AC_MSG_ERROR([could not find math library, are development libraries installed?])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], ,
    [AC_MSG_ERROR([could not find required version of GSL/CBLAS, are development libraries installed?])])
AC_CHECK_LIB([gsl], [gsl_spline_init], ,
    [AC_MSG_ERROR([could not find required version of GSL, are development libraries installed?])])

OPENMP_LDFLAGS="-fopenmp -lgomp"
LDFLAGS=${LDFLAGS}\ ${OPENMP_LDFLAGS}

dnl check for hdf5
AC_CHECK_HEADERS([hdf5.h], , [AC_MSG_ERROR([hdf5 header files not found, try CFLAGS=-I/path/to/headers])])
AC_CHECK_LIB([hdf5], [H5Acreate2], , [AC_MSG_ERROR([hdf5 library 1.8 not found, try LDFLAGS=-L/path/to/library])])

dnl check for test-runner, check
AM_CONDITIONAL([HAVE_CHECK], false)
PKG_CHECK_MODULES([check], [check >= 0.9.4], AM_CONDITIONAL([HAVE_CHECK], true), AC_MSG_WARN([check not found; testsuite not available]))

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h limits.h stdlib.h string.h sys/timeb.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AM_PROG_CC_C_O

# Checks for library functions.
AC_CHECK_FUNCS([floor memset pow sqrt strstr])

AC_CONFIG_FILES([Makefile
                 libLanlGeoMag/Lgm/Makefile
                 libLanlGeoMag/EopData/Makefile
                 libLanlGeoMag/Makefile
                 tests/Makefile
                 libLanlGeoMag/EopData/GetEopFiles
                 libLanlGeoMag/lgm.pc
                 lanlgeomag.spec
		 Tools/Makefile
                 ])
AC_OUTPUT
