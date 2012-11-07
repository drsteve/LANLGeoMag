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

AM_CPPFLAGS=-DDATADIR='$(datadir)'

DX_DOXYGEN_FEATURE(OFF)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(LanlGeoMag, Doxyfile, docs)


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

AC_ARG_ENABLE([static-tools], AS_HELP_STRING([--enable-static-tools], [Statically link tools (default: dynamic link). Implies --enable-tools]))
AS_IF([test "x$enable_static_tools" == "xyes"],
    [AM_CONDITIONAL([ENABLE_STATIC_TOOLS], true)
     AM_CONDITIONAL([ENABLE_TOOLS], true)
     enable_tools=yes
    ],
    [AM_CONDITIONAL([ENABLE_STATIC_TOOLS], false)]
)

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
dnl Need to directly specify SPICE location for MagEphem
dnl Can change this later to go hunting for them, but standard ld search path won't find it.
AC_ARG_WITH([cspice],
  [AS_HELP_STRING([--with-cspice=DIR], [Path to NAIF/cSPICE (contains include, lib).])],
  [cspice_CFLAGS=-I${withval}/include; cspice_LIBS=${withval}/lib/cspice.a] )
AC_ARG_WITH([cspice-includedir],
  [AS_HELP_STRING([--with-cspice-includedir=DIR], [Path to NAIF/cSPICE includes (contains SpiceUsr.h). Required if installing tools.])],
  [cspice_CFLAGS="-I$withval"]
)
AC_ARG_WITH([cspice-libdir],
  [AS_HELP_STRING([--with-cspice-libdir=DIR], [Path to NAIF/cSPICE library (contains cspice.a). Required if installing tools.])],
  [cspice_LIBS="${withval}/cspice.a"]
)
if test "x$enable_tools" == "xyes" && { test "x$cspice_CFLAGS" == "x" || test "x$cspice_LIBS" == "x"; }; then
   AC_MSG_FAILURE([--enable-tools requires either --with-cspice or both --with-cspice-libdir and --with-cspice-includedir])
fi
AC_SUBST(cspice_CFLAGS)
AC_SUBST(cspice_LIBS)

# Allow disabling of Qin-Denton files
AC_ARG_ENABLE([indices],
  [AS_HELP_STRING([--disable-indices], [Disable installation of geomagnetic index data (default: enable)])],
  [AM_CONDITIONAL([ENABLE_INDICES], [test x"$enableval" != x"no"])],
  [AM_CONDITIONAL([ENABLE_INDICES], true)]
)

# Check for python and SpacePy
build_py='no'
AC_ARG_WITH([python],
    [AS_HELP_STRING([--with-python=PYTHON],\
        [Name of (or path to) Python interpreter, or specify
         --without-python to skip Python bindings.])],
    [],
    [with_python=default]
)
if test x"${with_python}" == x"no"; then
    AC_MSG_NOTICE([skipping Python bindings (--without-python specified)])
else
    if test x"${with_python}" != x"default" && test x"${with_python}" != x"yes"; then
        PYTHON=${with_python}
    fi
    AM_PATH_PYTHON([2.6], , [:])
    dnl didn't find python, error or warn as appropriate
    if test x"$PYTHON" == x":"; then
        if test x"${with_python}" == x"default"; then
          AC_MSG_WARN(
          [Python 2.6+ not found; python support will not be built])
        else
            AC_MSG_FAILURE(
            [Python 2.6+ not found; cannot enable python support])
        fi
    else
	AC_MSG_CHECKING([for spacepy])
        ${PYTHON} -c "import spacepy" > /dev/null 2>/dev/null
        if test $? != 0; then
	    AC_MSG_RESULT([not found])
            if test x"${with_python}" == x"default"; then
                AC_MSG_WARN(
                  [spacepy not found; python support will not be built])
            else
                AC_MSG_FAILURE(
                [spacepy not found; cannot enable python support])
            fi
        else
	    AC_MSG_RESULT([found])
	    build_py="yes"
        fi
    fi
fi
AM_CONDITIONAL([HAVE_PYTHON], [test "$build_py" == "yes"])

#allow override of the Python output directory
AC_ARG_WITH([pythondir],
    AS_HELP_STRING([--with-pythondir=PATH],
    [Override the path to install Python modules [LIBDIR/pythonVER/site-packages/]]),
    [pythondir="$withval"], [])
AC_SUBST(pythondir)

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
AC_CHECK_LIB([perl], [perl_alloc], [AC_DEFINE([HAVE_LIBPERL], [1], [libperl found])], [AC_MSG_FAILURE([No libperl found, are perl development files installed?])])
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

# Check for openmp support
dnl note the linker does not need an explicit -lgomp, -fopenmp handles it
m4_ifdef([AC_OPENMP], [AC_OPENMP], [OPENMP_CFLAGS="-fopenmp"])
CFLAGS=${CFLAGS}\ ${OPENMP_CFLAGS}
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [return _OPENMP;])], , [AC_MSG_FAILURE([OpenMP support not found])])

# check for hdf5
AC_CHECK_HEADERS([hdf5.h], , [AC_MSG_FAILURE([hdf5 header files not found, try CFLAGS=-I/path/to/headers])])
AC_CHECK_LIB([hdf5], [H5Acreate2], , [AC_MSG_FAILURE([hdf5 library 1.8 not found, try LDFLAGS=-L/path/to/library])])

# check for argp (mostly for OS X)
AC_SEARCH_LIBS([argp_parse], [argp], , [AC_MSG_FAILURE([cannot find argp library])])

# check for test-runner, check
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
		 Python/Makefile
                 ])
AC_OUTPUT
