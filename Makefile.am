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

if DX_COND_pdf
    MAYBE_PDF = @DX_DOCDIR@/@DX_PROJECT@.pdf
    uninst_pdf = $(pdfdir)
endif
if DX_COND_ps
    MAYBE_PS = @DX_DOCDIR@/@DX_PROJECT@.ps
    uninst_ps = $(pdfdir)
endif
if DX_COND_html
    MAYBE_HTML = @DX_DOCDIR@/html
    uninst_html = $(htmldir)
endif
dist_pdf_DATA = $(MAYBE_PDF)
dist_ps_DATA = $(MAYBE_PS)
if DX_COND_doc
dist_noinst_DATA = @DX_DOCDIR@/@DX_PROJECT@.tag
endif

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
	$(MKDIR_P) -m 755 $(htmldir)/html/
	$(INSTALL_DATA) @DX_DOCDIR@/html/* $(htmldir)/html/

uninstall-local:
	-rm -rf $(htmldir)/html/
endif

if !DX_COND_doc
#if didn't ask for doxygen, just install whatever is already built
#(and uninstall whatever might be there)
uninst_html = $(htmldir)
uninst_ps = $(psdir)
uninst_pdf = $(pdfdir)

install-data-local:
	if [ -f @DX_DOCDIR@/html/index.html ]; then \
	    $(MKDIR_P) -m 755 $(htmldir)/html/; \
	    $(INSTALL_DATA) @DX_DOCDIR@/html/* $(htmldir)/html/; \
	fi
	if [ -f @DX_DOCDIR@/@DX_PROJECT@.ps ]; then \
	    $(MKDIR_P) -m 755 $(psdir); \
	    $(INSTALL_DATA) @DX_PROJECT@.ps $(psdir); \
	fi
	if [ -f @DX_DOCDIR@/@DX_PROJECT@.pdf ]; then \
	    $(MKDIR_P) -m 755 $(pdfdir); \
	    $(INSTALL_DATA) @DX_PROJECT@.pdf $(pdfdir); \
	fi

uninstall-local:
	-rm -rf $(htmldir)/html/
	-rm -f $(psdir)/@DX_PROJECT@.ps
	-rm -f $(pdfdir)/@DX_PROJECT@.pdf
endif

#remove any doc directories we've made
uninstall-hook:
	-rmdir $(uninst_ps) $(uninst_pdf) $(uninst_html) $(docdir)
