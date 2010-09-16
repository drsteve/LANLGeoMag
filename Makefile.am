#include $(top_srcdir)/aminclude.am

if HAVE_CHECK
SUBDIRS = libLanlGeoMag . tests
else
SUBDIRS = libLanlGeoMag
endif


EXTRA_DIST = Examples Doc ViewDriftShell

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -DDATADIR='"$(datadir)"'
