#make sure EopData is processed first, so on uninstall the data are gone.
SUBDIRS = Lgm EopData .

# Directory where Eop-related stuff will be installed
LGM_EOP_DATA_DIR   = $(datarootdir)/LanlGeoMag/EopData

# Directory where TS07-related files will be installed
LGM_TS07_DATA_DIR = $(datarootdir)/LanlGeoMag/Data/TS07D_FILES

# Directory where Index-related stuff will be installed
LGM_INDEX_DATA_DIR = $(datarootdir)/LanlGeoMag/Data

# Directory where JPL DE-related stuff will be installed
LGM_DE_DATA_DIR = $(datarootdir)/LanlGeoMag/Data

# Directory where Blue Marble 5400x2700 images will be installed
LGM_BLUEMARBLE_DATA_DIR = $(datarootdir)/LanlGeoMag/Images/BlueMarble

# Directory where other images (for ViewDriftSehll) will be installed
LGM_IMAGE_DATA_DIR = $(datarootdir)/LanlGeoMag/Images

# Directory where other images (for ViewDriftSehll) will be installed
LGM_STARS_DATA_DIR = $(datarootdir)/LanlGeoMag/Stars

pkgconfigdir       = $(prefix)/lib/pkgconfig
pkgconfig_DATA     = lgm.pc


if LINUX
    CRONDIR = /etc/cron.daily
else
    CRONDIR = /etc/periodic/daily
endif


#libdir                   = @prefix@/lib
lib_LTLIBRARIES          = libLanlGeoMag.la
libLanlGeoMag_la_SOURCES =  Lgm_AlphaOfK.c Lgm_DFI_RBF.c Lgm_Vec_RBF.c Lgm_B_FromScatteredData.c ComputeLstar.c DriftShell.c IntegralInvariant.c LFromIBmM.c \
	                        Lgm_B_internal.c Lgm_CTrans.c Lgm_DateAndTime.c Lgm_Eop.c Lgm_IGRF.c Lgm_InitMagInfo.c \
                            Lgm_MaxwellJuttner.c Lgm_Nutation.c Lgm_Octree.c Lgm_Quat.c Lgm_Sgp.c Lgm_SimplifiedMead.c  Lgm_SunPosition.c \
                            Lgm_Trace.c Lgm_TraceToEarth.c Lgm_TraceToSphericalEarth.c Lgm_Vec.c MagStep.c Lgm_QuadPack3.c \
                            Lgm_QuadPack.c Lgm_Cgm.c quicksort.c SbIntegral.c T87.c T89.c T89c.c TraceLine.c Lgm_TraceToMinBSurf.c  \
                            TraceToMinRdotB.c Lgm_TraceToMirrorPoint.c TraceToSMEquat.c T01S.c Tsyg_T01s.c T02.c Tsyg_T02.c TS04.c Tsyg2004.c \
                            Lgm_PolyRoots.c Lgm_SummersDiffCoeff.c Lgm_B_Dungey.c Tsyg2007.c TS07.c Tsyg1996.c T96.c TU82.c\
                            W.c Lgm_InitMagEphemInfo.c Lgm_AE8_AP8.c OP77.c OP88.c OlsenPfitzerDynamic.c OlsenPfitzerStatic.c IsoTimeStringToDateTime.c \
			                size.c Lgm_FluxToPsd.c xvgifwr2.c praxis.c Lgm_SphHarm.c Lgm_McIlwain_L.c Lgm_ElapsedTime.c Lgm_KdTree.c\
			                Lgm_ComputeLstarVersusPA.c Lgm_MagEphemWrite.c Lgm_MagEphemWriteHdf.c brent.c Lgm_CdipMirrorLat.c ComputeI_FromMltMlat.c ComputeI_FromMltMlat2.c \
                            Lgm_QinDenton.c Lgm_DiffCoeff_param.c Lgm_AE_index.c Lgm_Misc.c Lgm_HDF5.c Lgm_GradB.c Lgm_VelStep.c Lgm_Utils.c DynamicMemory.h \
			                Lgm_Metadata.c  Lgm_PriorityQueue.c TraceToYZPlane.c Lgm_InitNrlMsise00.c Lgm_NrlMsise00.c Lgm_Coulomb.c\
			                Lgm_Ellipsoid.c Lgm_DipEquator.c \
                            Lgm_JPLephem.c  Lgm_Eclipse.c Lgm_TabularBessel.c



libLanlGeoMag_la_LDFLAGS = $(AM_LDFLAGS) @PERL_LDFLAGS@ @OPENMP_CFLAGS@ 
libLanlGeoMag_la_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ @OPENMP_CFLAGS@ -Wno-incompatible-pointer-types
libLanlGeoMag_la_CPPFLAGS = $(AM_CPPFLAGS) -DDATADIR=\"$(datadir)\" -DHASH_FUNCTION=HASH_SFH -DHASH_BLOOM=32 @PERL_CPPFLAGS@

#EXTRA_DIST = TS07D_FILES DE_FILES Images Stars
EXTRA_DIST = TS07D_FILES DE_FILES Stars

install-data-hook:
if ENABLE_INDICES
	$(MKDIR_P) $(DESTDIR)$(LGM_TS07_DATA_DIR)
	$(RSYNC_CMD) -arv $(srcdir)/TS07D_FILES/* $(DESTDIR)$(LGM_TS07_DATA_DIR)/.
	$(MKDIR_P) $(DESTDIR)$(LGM_INDEX_DATA_DIR)
	@echo Qin-Denton index files should be placed in $(DESTDIR)$(LGM_INDEX_DATA_DIR)
	$(MKDIR_P) $(DESTDIR)$(LGM_DE_DATA_DIR)
	$(RSYNC_CMD) -arv $(srcdir)/DE_FILES $(DESTDIR)$(LGM_DE_DATA_DIR)/.
	@echo Additional DE files should be placed in $(DESTDIR)$(LGM_DE_DATA_DIR)
	$(MKDIR_P) $(DESTDIR)$(LGM_STARS_DATA_DIR)
	$(RSYNC_CMD) -arv $(srcdir)/Stars/* $(DESTDIR)$(LGM_STARS_DATA_DIR)/.
endif
#	# $(MKDIR_P) $(DESTDIR)$(LGM_IMAGE_DATA_DIR)
#	# $(RSYNC_CMD) -arv $(srcdir)/Images/* $(DESTDIR)$(LGM_IMAGE_DATA_DIR)/. # should cover BlueMarble stuff also.

uninstall-hook:
if ENABLE_INDICES
endif
	-rmdir $(DESTDIR)$(LGM_INDEX_DATA_DIR) #don't delete user data
	-chmod -R 777 $(DESTDIR)$(LGM_TS07_DATA_DIR)
	-rm -rf $(DESTDIR)$(LGM_TS07_DATA_DIR)
	-chmod -R 777 $(DESTDIR)$(LGM_DE_DATA_DIR)
	-rm -rf $(DESTDIR)$(LGM_DE_DATA_DIR)
	-chmod -R 777 $(DESTDIR)$(LGM_IMAGE_DATA_DIR)
	-rm -rf $(DESTDIR)$(LGM_IMAGE_DATA_DIR)
	-chmod -R 777 $(DESTDIR)$(LGM_STARS_DATA_DIR)
	-rm -rf $(DESTDIR)$(LGM_STARS_DATA_DIR)
	-rmdir $(DESTDIR)$(datarootdir)/LanlGeoMag
    # 	echo "WILDCARD: "
    # 	echo $(wildcard \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/tailamhr_?_??.par \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/tailamebhr?.par \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/TAIL_PAR.zip \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/Coeffs/????_???/2008_???_??_??.par \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/DE_FILES/jpl_de421.h5 \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Stars/Stars.txt)
    # 	for fname in $(wildcard \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/tailamhr_?_??.par \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/tailamebhr?.par \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/TAIL_PAR/TAIL_PAR.zip \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/TS07D_FILES/Coeffs/????_???/2008_???_??_??.par \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Data/DE_FILES/jpl_de421.h5 \
    # $(DESTDIR)$(prefix)/share/LanlGeoMag/Stars/Stars.txt); \
    # 	do \
    # 	    echo "fname: $${fname}"; \
    # 	    chmod 644 $${fname}; \
    # 	    /bin/rm -f $${fname}; \
    # 	done
