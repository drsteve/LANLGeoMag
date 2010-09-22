#include $(top_srcdir)/aminclude.am

if HAVE_CHECK
SUBDIRS = libLanlGeoMag . tests
else
SUBDIRS = libLanlGeoMag
endif

ALL_EXAMPLE_FILES = $(shell cat example_file_list)
ALL_DOC_FILES = Doc/GSM.png Doc/CoordFig.eps Doc/GSM_SM_FIG.svg \
	Doc/Geodetic.png Doc/CoordFig2.svg Doc/CoordFig.svg Doc/GSM_SM_FIG.png \
	Doc/SM.png Doc/SunPosition.png Doc/SunPosition.eps \
	Doc/LanlGeoMag_Documentation.odt Doc/LanlGeoMag_Documentation.pdf \
	Doc/CoordFig2.png Doc/IridiumFlareFrame.png Doc/DUT1_DAT.png \
	Doc/Geodetic.svg Doc/GSM_SM_FIG.eps.gz Doc/CoordFig.png \
	Doc/EopParams_1990_2010.png Doc/TimeDiffs.png
ALL_VIEW_DRIFT_SHELL_FILES = $(shell cat view_drift_shell_file_list)

#EXTRA_DIST = Examples Doc ViewDriftShell
EXTRA_DIST = $(ALL_EXAMPLE_FILES) $(ALL_DOC_FILES) \
	$(ALL_VIEW_DRIFT_SHELL_FILES) example_file_list \
	view_drift_shell_file_list

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -DDATADIR='"$(datadir)"'
