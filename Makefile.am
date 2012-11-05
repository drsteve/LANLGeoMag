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


include aminclude.am                                                                                                                                                                                             

if DX_COND_man                                                                                                                                                                                                   

# You'd probably want to post-process man pages and installed the patched                                                                                                                                        
# versions.                                                                                                                                                                                                      

#dist_man1_MANS = @DX_DOCDIR@/man/man1/doxample.c.1                                                                                                                                                               
$(dist_man1_MANS): doxygen-doc                                                                                                                                                                                   
#$(man3_MANS): doxygen-doc

endif                                           

MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = lanlgeomag.spec $(ALL_EXAMPLE_FILES) $(ALL_DOC_FILES) \
	Python HACKING $(DX_CONFIG)

ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -DDATADIR='"$(datadir)"'
