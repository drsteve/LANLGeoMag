#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([LanlGeoMag], [1.5.16], [mghenderson@lanl.gov])
AC_CONFIG_SRCDIR([libLanlGeoMag/quicksort.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([tar-ustar])


# Checks for programs.
AC_PROG_CC
AC_PROG_SED
AC_PROG_MKDIR_P

AC_PROG_LIBTOOL
#LT_INIT

PKG_INSTALLDIR

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h limits.h stdlib.h string.h sys/timeb.h unistd.h math.h])



AC_PREREQ([2.59])
THE_PKG_RELEASE=8dio dnl this affects the RPM package release #
PKG_RELEASE=$THE_PKG_RELEASE
AC_SUBST([PKG_RELEASE])

AC_LANG([C])
#AC_DEFINE([LGM_EOP_DATA_DIR],["/usr/local/share/LanlGeoMag/EopData"],["Define directory where Eop data can be found"])
#AH_VERBATIM([PDDH], [#define LGM_EOP_DATA_DIR "@datadir@/LanlGeoMag/EopData"])
#AC_DEFINE([PDD], ["PDD @prefix@/data"], ["dude"])
AC_DEFINE([LGM_EOP_DATA_DIR],[(DATADIR "/LanlGeoMag/EopData")],[Define directory where Eop data can be found])
AC_DEFINE([LGM_INDEX_DATA_DIR],[(DATADIR "/LanlGeoMag/Data")],[Define directory where Geophysical Index Data (e.g. Kp, AE, Dst, etc. data can be found])
AC_DEFINE([LGM_TS07_DATA_DIR],[(DATADIR "/LanlGeoMag/Data/TS07D_FILES")],[Define directory where TS07 data files are located])
#AM_INIT_AUTOMAKE



AM_CPPFLAGS=-DDATADIR='$(datadir)'

DX_DOXYGEN_FEATURE(OFF)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(LanlGeoMag, Doxyfile, Doc)

dnl we need an overall test of whether doxyen was enabled in the first
dnl place
AC_DEFUN([DX_TEST_OVERALL], [test "DX_DOXYGEN_FEATURE" = 1])
AM_CONDITIONAL(DX_OVERALL, DX_TEST_OVERALL])

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


#check for rsync (but just use cp if don't have it)
AC_PATH_PROG([RSYNC_CMD], [rsync], [cp])

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

if test "x$enable_tools" == "xyes"; then
   OLD_CPPFLAGS=${CPPFLAGS}
   CPPFLAGS=${CPPFLAGS}\ ${cspice_CFLAGS}
   OLD_LIBS=${LIBS}
   LIBS=${LIBS}\ ${cspice_LIBS}
   AC_CHECK_HEADER([SpiceUsr.h], , [AC_MSG_FAILURE([SpiceUsr.h not found or not usable.])])
   AC_MSG_CHECKING([for cspice.a])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [return 0;])], [AC_MSG_RESULT([found])], [AC_MSG_FAILURE([cspice.a not found or not usable.])])
   CPPFLAGS=${OLD_CPPFLAGS}
   LIBS=${OLD_LIBS}
fi

# Allow disabling of TS07D files
AC_ARG_ENABLE([indices],
  [AS_HELP_STRING([--disable-indices], [Disable installation of geomagnetic index data for TS07D (default: enable)])],
  [AM_CONDITIONAL([ENABLE_INDICES], [test x"$enableval" != x"no"])],
  [AM_CONDITIONAL([ENABLE_INDICES], true)]
)

# Allow disabling of ALL documentation
AC_ARG_ENABLE([docs],
  [AS_HELP_STRING([--disable-docs], [Disable installation of pre-built documentation. (doxygen docs will still be installed if --enable-doxygen-docs is specified.)])],
  [AM_CONDITIONAL([ENABLE_DOCS], [test x"$enableval" != x"no"])],
  [AM_CONDITIONAL([ENABLE_DOCS], true)]
)

# Check if want call graphs
DX_ENABLE_CALL_GRAPH=YES
AC_ARG_ENABLE([doxygen-call-graph],
  [AS_HELP_STRING([--disable-doxygen-call-graph], [Disable call graph if building doxygen docs.])],
  [AS_IF([test x"$enableval" == x"no"], [DX_ENABLE_CALL_GRAPH=NO])]
)
DX_ENV_APPEND(CALL_GRAPH, $DX_ENABLE_CALL_GRAPH)

# Check if want PDF linkx
DX_ENABLE_PDF_HYPERLINKS=YES
AC_ARG_ENABLE([doxygen-pdf-links],
  [AS_HELP_STRING([--disable-doxygen-pdf-links], [Disable links in PDF if building doxygen docs. May help with strange LaTeX errors])],
  [AS_IF([test x"$enableval" == x"no"], [DX_ENABLE_PDF_HYPERLINKS=NO])]
)
DX_ENV_APPEND(PDF_LINKS, $DX_ENABLE_PDF_HYPERLINKS)


# Check for check
#AM_PATH_CHECK(0.9.4,[have_check="yes"], AC_MSG_WARN([Check not found; cannot run unit tests!]) [have_check="no"])
#AM_CONDITIONAL(HAVE_CHECK, test x"$have_check" = "xyes")





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

# Make PERL dependency optional
AC_ARG_ENABLE([perl], [AS_HELP_STRING([--enable-perl],\
    [Enables additional PERL functionality in LanlGeoMag])])
AS_IF([test "x$enableval" != "xyes"],
    [ENABLE_PERL=0],
    [ENABLE_PERL=1]
)

# Check for PERL if --with-perl specified
if test x"${ENABLE_PERL}" == x"1"; then
   AC_CHECK_PROG([HAVE_PERL], [perl], true, false)
   if test x"${HAVE_PERL}" == x"false"; then
       AC_MSG_FAILURE([Unable to run perl])
   fi
fi

# Define the PERL flags if using PERL
if test x"${HAVE_PERL}" == x"true"; then
   AX_PROG_PERL_MODULES( [ExtUtils::Embed], , [AC_MSG_FAILURE(
		      	 Required Perl module ExtUtils::Embed not found)])
   dnl Get flags for embedded perl and make sure it works
   dnl Save old CFLAGS/LDFLAGS so can restore after the check
   PERL_CFLAGS=`perl -MExtUtils::Embed -e ccopts`
   PERL_LDFLAGS=`perl -MExtUtils::Embed -e ldopts`
   #get out JUST the rpath bit, since this needs to go in LDFLAGS
   PERL_RPATH=`echo ${PERL_LDFLAGS} | ${SED} 's;.*\(\-Wl,\-rpath[[^ ]]*\).*;\1;'`
fi

HDF5_LDFLAGS=`pkg-config --libs hdf5 2>/dev/null`
HDF5_CFLAGS=`pkg-config --cflags hdf5 2>/dev/null`

CHECK_LIBS=`pkg-config --libs check 2>/dev/null`
CHECK_CFLAGS=`pkg-config --cflags check 2>/dev/null`

OLD_CFLAGS=${CFLAGS}
OLD_LDFLAGS=${LDFLAGS}
OLD_LIBS=${LIBS}
CFLAGS=${CFLAGS}\ ${PERL_CFLAGS}\ ${HDF5_CFLAGS}
LDFLAGS=${LDFLAGS}\ ${PERL_LDFLAGS}\ ${HDF5_LDFLAGS}

# Check for the perl library after setting PERL_LDFLAGS
if test x"${HAVE_PERL}" == x"true"; then
   AC_CHECK_LIB([perl], [perl_alloc],\
   [AC_DEFINE([HAVE_LIBPERL], [1], [libperl found])],\
   [AC_MSG_FAILURE([No libperl found, are perl development files installed?])])
fi

CFLAGS=${OLD_CFLAGS}
LDFLAGS=${OLD_LDFLAGS}\ ${PERL_RPATH}
LIBS=${OLD_LIBS}
PERL_CPPFLAGS=-DENABLE_PERL=${ENABLE_PERL}
AC_SUBST(PERL_CPPFLAGS)
AC_SUBST(PERL_CFLAGS)
AC_SUBST(PERL_LDFLAGS)
AC_SUBST(HDF5_CFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(CHECK_CFLAGS)
AC_SUBST(CHECK_LIBS)

dnl Check for gsl, http://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html
AC_CHECK_LIB([m], [cos], ,
    [AC_MSG_ERROR([could not find math library, are development libraries installed?])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], ,
    [AC_MSG_ERROR([could not find required version of GSL/CBLAS, are development libraries installed?])])
AC_CHECK_LIB([gsl], [gsl_spline_init], ,
    [AC_MSG_ERROR([could not find required version of GSL, are development libraries installed?])])

# Check for openmp support
m4_ifdef([AC_OPENMP], [
dnl use the nice macro if we have it
AC_OPENMP
AS_IF([test "x$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" == "xunsupported"], [
    USE_OPENMP=0], [
    USE_OPENMP=1])
], [
dnl don't have the nice macro, just try one quick flag...
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
OLD_CFLAGS=${CFLAGS}
CFLAGS=${CFLAGS}\ -fopenmp -lgomp
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [return _OPENMP;])], [
USE_OPENMP=1
OPENMP_CFLAGS="-fopenmp -lgomp"
], [
USE_OPENMP=0])
CFLAGS=${OLD_CFLAGS}])
# if we don't have openmp support, check if we asked for it...
if test $USE_OPENMP -eq 0; then
    if test "x$enable_openmp" = "xyes"; then
        AC_MSG_FAILURE([OpenMP support not found.])
    fi
    if test "x$enable_openmp" != "xno"; then
        AC_MSG_WARN([OpenMP support not found; building without OpenMP])
    fi
fi
AC_SUBST([OPENMP_CFLAGS])
AC_DEFINE_UNQUOTED([USE_OPENMP], [$USE_OPENMP], [
Enable multithreading/processing w/OpenMP])

# check for hdf5
AM_CONDITIONAL([HAVE_HDF5], false)
PKG_CHECK_MODULES([hdf5], [hdf5], AM_CONDITIONAL([HAVE_HDF5], true), [
AC_MSG_WARN([hdf5 not found in pkgconfig, hopefully it's in normal paths])])
#Can't do the assignment until out of the PKG_CHECK
AM_COND_IF(HAVE_HDF5, [
LIBS=${LIBS}\ $hdf5_LIBS
CFLAGS=${CFLAGS}\ $hdf5_CFLAGS])

AC_CHECK_HEADERS([hdf5.h], , [AC_MSG_FAILURE([hdf5 header files not found, try CFLAGS=-I/path/to/headers])])
AC_CHECK_LIB([hdf5], [H5Acreate2], , [AC_MSG_FAILURE([hdf5 library 1.8 not found, try LDFLAGS=-L/path/to/library])])

# check for argp (mostly for OS X)
AC_SEARCH_LIBS([argp_parse], [argp], , [AC_MSG_FAILURE([cannot find argp library])])

# check for test-runner, check
AM_CONDITIONAL([HAVE_CHECK], false)
PKG_CHECK_MODULES([check], [check >= 0.9.4], AM_CONDITIONAL([HAVE_CHECK], true), [
OLD_LIBS=${LIBS}
LIBS=${LIBS}\ -lcheck
AC_MSG_CHECKING([for alternate check])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <check.h>], [])], [
AC_MSG_RESULT([found])
AM_CONDITIONAL([HAVE_CHECK], true)], [
AC_MSG_RESULT([not found])
AC_MSG_WARN([check not found; testsuite not available])
LIBS=${OLD_LIBS}])])


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
                 Doc/LanlGeoMag.tag
                 ])
AC_OUTPUT
