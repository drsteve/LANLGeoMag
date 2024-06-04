if HAVE_CHECK
    MAYBE_TESTS = tests
endif
if ENABLE_TOOLS
    MAYBE_TOOLS = Tools
endif
if HAVE_PYTHON
    MAYBE_PY = Python
endif
#SUBDIRS = libLanlGeoMag . $(MAYBE_TOOLS) $(MAYBE_TESTS) $(MAYBE_PY)
SUBDIRS = libLanlGeoMag . $(MAYBE_TOOLS) $(MAYBE_TESTS) Python

if DX_OVERALL
    if DX_COND_pdf
        MAYBE_PDF = @DX_DOCDIR@/@DX_PROJECT@.pdf
        uninst_pdf = $(pdfdir)
    endif
endif
if DX_COND_ps
    MAYBE_PS = @DX_DOCDIR@/@DX_PROJECT@.ps
    uninst_ps = $(psdir)
endif
if DX_COND_html
    MAYBE_HTML = @DX_DOCDIR@/html
    uninst_html = $(htmldir)
endif
if ENABLE_DOCS
    MAYBE_STATIC_PDF = Doc/LanlGeoMag_Documentation.pdf
endif
if DX_COND_pdf
    dist_pdf_DATA = $(MAYBE_STATIC_PDF)
else
    dist_pdf_DATA = $(MAYBE_STATIC_PDF)
endif

nodist_ps_DATA = $(MAYBE_PS)
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

DEBIAN_EXTRA_DIST = debian/rules \
    debian/lanlgeomag.install \
    debian/copyright \
    debian/compat \
    debian/changelog \
    debian/control \
    debian/lanlgeomag-dev.dirs \
    debian/source/format \
    debian/files \
    debian/lanlgeomag.dirs \
    debian/lanlgeomag-dev.install

EXTRA_DIST = README.md lanlgeomag.spec $(ALL_EXAMPLE_FILES) $(ALL_DOC_FILES) \
	Python HACKING $(DX_CONFIG) $(MAYBE_HTML) fix_doxygen.sed  $(DEBIAN_EXTRA_DIST)

ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -DDATADIR='"$(datadir)"' -Wno-incompatible-pointer-types

README: README.md
	pandoc -f markdown -t plain --wrap=none $< -o $@

CLEANFILES = README

if DX_COND_html
#and make sure html is build when doing a make all
#(for some reason this bit doesn't honor DESTDIR without doing it explicitly)
all-local: @DX_DOCDIR@/html


uninstall-local:
	-rm -rf $(htmldir)/html/
endif

if !DX_COND_doc
if ENABLE_DOCS
#if didn't ask for doxygen, but want docs, just install whatever is already built
#(and uninstall whatever might be there)
uninst_html = $(htmldir)
uninst_ps = $(psdir)
uninst_pdf = $(pdfdir)

install-data-local:
if DX_OVERALL
if DX_COND_html
	$(MKDIR_P) -m 755 $(DESTDIR)$(htmldir)/html/
	$(INSTALL_DATA) @DX_DOCDIR@/html/*.* $(DESTDIR)$(htmldir)/html/
	if [ -d @DX_DOCDIR@/html/search ]; then \
	    $(MKDIR_P) -m 755 $(DESTDIR)$(htmldir)/html/search; \
	    $(INSTALL_DATA) @DX_DOCDIR@/html/search/*.* $(DESTDIR)$(htmldir)/html/search; \
	fi
	if [ -f @DX_DOCDIR@/html/index.html ]; then \
	    $(MKDIR_P) -m 755 $(htmldir)/html/; \
	    $(INSTALL_DATA) @DX_DOCDIR@/html/*.* $(htmldir)/html/; \
	fi
	if [ -d @DX_DOCDIR@/html/search ]; then \
	    $(MKDIR_P) -m 755 $(htmldir)/html/search; \
	    $(INSTALL_DATA) @DX_DOCDIR@/html/search/*.* $(htmldir)/html/search; \
	fi
endif				# DX_COND_html
	if [ -f @DX_DOCDIR@/@DX_PROJECT@.ps ]; then \
	    $(MKDIR_P) -m 755 $(psdir); \
	    $(INSTALL_DATA) @DX_DOCDIR@/@DX_PROJECT@.ps $(psdir); \
	fi
	if [ -f $(MAYBE_PDF) ]; then \
	    $(MKDIR_P) -m 755 $(pdfdir); \
	    $(INSTALL_DATA) $(MAYBE_PDF) $(pdfdir); \
	fi
endif				# DX_OVERALL
uninstall-local:
	-rm -rf $(htmldir)/html/
	-rm -f $(psdir)/@DX_PROJECT@.ps
	-rm -f $(pdfdir)/`basename $(MAYBE_PDF)`
endif ENABLE_DOCS
endif !DX_COND_doc

#remove any doc directories we've made
uninstall-hook:
	-rmdir $(uninst_ps) $(uninst_pdf) $(uninst_html) $(docdir)

# CLEANFILES += \
#     ./share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/tailamhr_?_??.par \
#     ./share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/tailamebhr?.par \
#     ./share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/TAIL_PAR.zip \
#     ./share/LanlGeoMag/Data/TS07D_FILES/Coeffs/????_???/2008_???_??_??.par \
#     ./share/LanlGeoMag/Data/DE_FILES/jpl_de421.h5 \
#     ./share/LanlGeoMag/Stars/Stars.txt
