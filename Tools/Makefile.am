LGMSRCDIR = $(top_srcdir)/libLanlGeoMag/
bin_PROGRAMS = MagEphemFromSpiceKernel
MagEphemFromSpiceKernel_SOURCES = MagEphemFromSpiceKernel.c
MagEphemFromSpiceKernel_LDADD = $(top_builddir)/libLanlGeoMag/libLanlGeoMag.la @cspice_LIBS@
if ENABLE_STATIC_TOOLS
    MagEphemFromSpiceKernel_LDFLAGS = $(AM_LDFLAGS) @PERL_LDFLAGS@ -static @OPENMP_CFLAGS@
    MagEphemFromSpiceKernel_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ -I$(LGMSRCDIR) -I$(LGMSRCDIR)/Lgm @cspice_CFLAGS@ @OPENMP_CFLAGS@
else
    MagEphemFromSpiceKernel_CFLAGS = $(AM_CFLAGS) -I$(LGMSRCDIR) -I$(LGMSRCDIR)/Lgm @cspice_CFLAGS@
endif

#CFLAGS       = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
#SPICE_LIBS = '/home/mgh/SPICE/cspice/lib/cspice.a'
#SPICE_INC  = '-I/home/mgh/SPICE/cspice/include'
#HDF5_LIBS  = '-L/packages/lib/hdf5-1.8.7'
#HDF5_INC   = '-I/packages/include/hdf5-1.8'
#LIBS       = `pkg-config --libs lgm`
#CFLAGS     = `pkg-config --cflags lgm` $(SPICE_INC)
