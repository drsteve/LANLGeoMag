if HAVE_CHECK
    MAYBE_TESTS = tests
endif
if ENABLE_TOOLS
    MAYBE_TOOLS = Tools
endif
if HAVE_PYTHON
    MAYBE_PY = Python
endif
SUBDIRS = libLanlGeoMag . $(MAYBE_TOOLS) $(MAYBE_TESTS) $(MAYBE_PY)

documentationdir = $(datarootdir)/LanlGeoMag/docs
if DX_COND_pdf
    MAYBE_PDF = @DX_DOCDIR@/@DX_PROJECT@.pdf
endif
if DX_COND_ps
    MAYBE_PS = @DX_DOCDIR@/@DX_PROJECT@.ps
endif
if DX_COND_html
    MAYBE_HTML = @DX_DOCDIR@/html
endif
dist_documentation_DATA = $(MAYBE_PDF) $(MAYBE_PS)
dist_noinst_DATA = @DX_DOCDIR@/@DX_PROJECT@.tag

GZIP_ENV = -7

ALL_EXAMPLE_FILES = Examples
ALL_DOC_FILES = Doc/GSM.png Doc/CoordFig.eps Doc/GSM_SM_FIG.svg \
	Doc/Geodetic.png Doc/CoordFig2.svg Doc/CoordFig.svg Doc/GSM_SM_FIG.png \
	Doc/SM.png Doc/SunPosition.png Doc/SunPosition.eps \
	Doc/LanlGeoMag_Documentation.odt Doc/LanlGeoMag_Documentation.pdf \
	Doc/CoordFig2.png Doc/IridiumFlareFrame.png Doc/DUT1_DAT.png \
	Doc/Geodetic.svg Doc/GSM_SM_FIG.eps.gz Doc/CoordFig.png \
	Doc/EopParams_1990_2010.png Doc/TimeDiffs.png
#ALL_VIEW_DRIFT_SHELL_FILES = ViewDriftShell


include doxygen_include.am

MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = lanlgeomag.spec $(ALL_EXAMPLE_FILES) $(ALL_DOC_FILES) \
	Python HACKING $(DX_CONFIG) $(MAYBE_HTML) fix_doxygen.sed

ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -DDATADIR='"$(datadir)"'

if DX_COND_html
#and make sure html is build when doing a make all
all-local: @DX_DOCDIR@/html

install-data-local:
	$(MKDIR_P) -m 755 $(documentationdir)/html/
	$(INSTALL_DATA) @DX_DOCDIR@/html/* $(documentationdir)/html/

uninstall-local:
	-rm -rf $(documentationdir)/html/

uninstall-hook:
	-rmdir $(documentationdir)
	-rmdir $(datarootdir)/LanlGeoMag/
endif
