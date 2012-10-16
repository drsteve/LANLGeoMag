LGMSRCDIR=../libLanlGeoMag/
bin_PROGRAMS = MagEphemFromSpiceKernel
MagEphemFromSpiceKernel_SOURCES = MagEphemFromSpiceKernel.c
MagEphemFromSpiceKernel_LDADD = $(LGMSRCDIR)/libLanlGeoMag.la /packages/lib/cspice.a #@cspice_LIBS@
#hard-code the include path, since the headers don't go anywhere funny pre-install.
#This may change if using, say, autoheader.
MagEphemFromSpiceKernel_CFLAGS = $(AM_CFLAGS) -I$(LGMSRCDIR) -I$(LGMSRCDIR)/Lgm -I/packages/include/naif #@cspice_CFLAGS@

#CFLAGS       = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
#SPICE_LIBS = '/home/mgh/SPICE/cspice/lib/cspice.a'
#SPICE_INC  = '-I/home/mgh/SPICE/cspice/include'
#HDF5_LIBS  = '-L/packages/lib/hdf5-1.8.7'
#HDF5_INC   = '-I/packages/include/hdf5-1.8'
#LIBS       = `pkg-config --libs lgm`
#CFLAGS     = `pkg-config --cflags lgm` $(SPICE_INC)
