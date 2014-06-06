LGMSRCDIR = $(top_srcdir)/libLanlGeoMag/

bin_PROGRAMS = MagEphemFromSpiceKernel MagEphemFromTLE LastClosedDriftShell

MagEphemFromSpiceKernel_SOURCES = MagEphemFromSpiceKernel.c
MagEphemFromSpiceKernel_LDADD = $(top_builddir)/libLanlGeoMag/libLanlGeoMag.la @cspice_LIBS@
if ENABLE_STATIC_TOOLS
    MagEphemFromSpiceKernel_LDFLAGS = $(AM_LDFLAGS) @PERL_LDFLAGS@ -static @OPENMP_CFLAGS@
    MagEphemFromSpiceKernel_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ @cspice_CFLAGS@ @OPENMP_CFLAGS@
else
    MagEphemFromSpiceKernel_LDFLAGS = $(AM_LDFLAGS) @OPENMP_CFLAGS@
    MagEphemFromSpiceKernel_CFLAGS = $(AM_CFLAGS) @cspice_CFLAGS@ @OPENMP_CFLAGS@
endif
#http://www.gnu.org/software/automake/manual/html_node/Flag-Variables-Ordering.html
MagEphemFromSpiceKernel_CPPFLAGS = -I$(LGMSRCDIR) -I$(LGMSRCDIR)/Lgm $(AM_CPPFLAGS)

MagEphemFromTLE_SOURCES = MagEphemFromTLE.c
MagEphemFromTLE_LDADD = $(top_builddir)/libLanlGeoMag/libLanlGeoMag.la @cspice_LIBS@
if ENABLE_STATIC_TOOLS
    MagEphemFromTLE_LDFLAGS = $(AM_LDFLAGS) @PERL_LDFLAGS@ -static @OPENMP_CFLAGS@
    MagEphemFromTLE_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ @cspice_CFLAGS@ @OPENMP_CFLAGS@
else
    MagEphemFromTLE_LDFLAGS = $(AM_LDFLAGS) @OPENMP_CFLAGS@
    MagEphemFromTLE_CFLAGS = $(AM_CFLAGS) @cspice_CFLAGS@ @OPENMP_CFLAGS@
endif
#http://www.gnu.org/software/automake/manual/html_node/Flag-Variables-Ordering.html
MagEphemFromTLE_CPPFLAGS = -I$(LGMSRCDIR) -I$(LGMSRCDIR)/Lgm $(AM_CPPFLAGS)

LastClosedDriftShell_SOURCES = LastClosedDriftShell.c
LastClosedDriftShell_LDADD = $(top_builddir)/libLanlGeoMag/libLanlGeoMag.la @cspice_LIBS@
if ENABLE_STATIC_TOOLS
    LastClosedDriftShell_LDFLAGS = $(AM_LDFLAGS) @PERL_LDFLAGS@ -static @OPENMP_CFLAGS@
    LastClosedDriftShell_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ @OPENMP_CFLAGS@
else
    LastClosedDriftShell_LDFLAGS = $(AM_LDFLAGS) @OPENMP_CFLAGS@
    LastClosedDriftShell_CFLAGS = $(AM_CFLAGS) @OPENMP_CFLAGS@
endif
#http://www.gnu.org/software/automake/manual/html_node/Flag-Variables-Ordering.html
LastClosedDriftShell_CPPFLAGS = -I$(LGMSRCDIR) -I$(LGMSRCDIR)/Lgm $(AM_CPPFLAGS)
