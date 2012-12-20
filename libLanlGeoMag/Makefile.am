#make sure EopData is processed first, so on uninstall the data are gone.
SUBDIRS = Lgm EopData .

# Directory where Eop-related stuff will be installed
LGM_EOP_DATA_DIR   = $(datarootdir)/LanlGeoMag/EopData

# Directory where TS07-related files will be installed
LGM_TS07_DATA_DIR = $(datarootdir)/LanlGeoMag/Data

# Directory where Index-related stuff will be installed
LGM_INDEX_DATA_DIR = $(datarootdir)/LanlGeoMag/Data

pkgconfigdir       = $(libdir)/pkgconfig
pkgconfig_DATA     = lgm.pc


if LINUX
    CRONDIR = /etc/cron.daily
else
    CRONDIR = /etc/periodic/daily
endif


#libdir                   = @prefix@/lib
lib_LTLIBRARIES          = libLanlGeoMag.la
libLanlGeoMag_la_SOURCES =  Lgm_AlphaOfK.c Lgm_DFI_RBF.c Lgm_B_FromScatteredData.c ComputeLstar.c DriftShell.c IntegralInvariant.c LFromIBmM.c \
	                        Lgm_B_internal.c Lgm_CTrans.c Lgm_DateAndTime.c Lgm_Eop.c Lgm_IGRF.c Lgm_InitMagInfo.c \
                            Lgm_MaxwellJuttner.c Lgm_Nutation.c Lgm_Octree.c Lgm_Quat.c Lgm_Sgp.c Lgm_SimplifiedMead.c  Lgm_SunPosition.c \
                            Lgm_Trace.c Lgm_TraceToEarth.c Lgm_TraceToSphericalEarth.c Lgm_Vec.c MagStep.c Lgm_QuadPack3.c \
                            Lgm_QuadPack.c Lgm_Cgm.c quicksort.c SbIntegral.c T87.c T89.c T89c.c TraceLine.c Lgm_TraceToMinBSurf.c  \
                            TraceToMinRdotB.c Lgm_TraceToMirrorPoint.c TraceToSMEquat.c T01S.c Tsyg_T01s.c TS04.c Tsyg2004.c \
                            Lgm_PolyRoots.c Lgm_SummersDiffCoeff.c Lgm_B_Dungey.c Tsyg2007.c TS07.c\
                            W.c Lgm_InitMagEphemInfo.c Lgm_AE8_AP8.c OP77.c OlsenPfitzerStatic.c IsoTimeStringToDateTime.c \
			                size.c Lgm_FluxToPsd.c xvgifwr2.c praxis.c Lgm_SphHarm.c Lgm_McIlwain_L.c Lgm_ElapsedTime.c \
			                Lgm_ComputeLstarVersusPA.c Lgm_MagEphemWrite.c Lgm_MagEphemWriteHdf.c brent.c Lgm_CdipMirrorLat.c ComputeI_FromMltMlat.c \
                            Lgm_QinDenton.c Lgm_DiffCoeff_param.c Lgm_AE_index.c Lgm_Misc.c Lgm_HDF5.c Lgm_GradB.c Lgm_VelStep.c Lgm_Utils.c DynamicMemory.h \
			    Lgm_Metadata.c  #Lgm_Enkf_1D_RadDiff.c #Lgm_FastPowPoly.c

libLanlGeoMag_la_LDFLAGS = $(AM_LDFLAGS) @OPENMP_CFLAGS@
libLanlGeoMag_la_LIBADD = @PERL_LDFLAGS@
libLanlGeoMag_la_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ @OPENMP_CFLAGS@
libLanlGeoMag_la_CPPFLAGS = $(AM_CPPFLAGS) -DDATADIR=\"$(datadir)\" -DHASH_FUNCTION=HASH_SFH -DHASH_BLOOM=32

EXTRA_DIST = QinDenton TS07D_FILES

install-data-hook:
if ENABLE_INDICES
	$(MKDIR_P) $(DESTDIR)$(LGM_INDEX_DATA_DIR)
	$(RSYNC_CMD) -a $(srcdir)/QinDenton $(DESTDIR)$(LGM_INDEX_DATA_DIR)/.
	$(RSYNC_CMD) -a $(srcdir)/TS07D_FILES $(DESTDIR)$(LGM_TS07_DATA_DIR)/.
endif

uninstall-hook:
if ENABLE_INDICES
	-rm -rf $(DESTDIR)$(LGM_INDEX_DATA_DIR)
	-rm -rf $(DESTDIR)$(LGM_TS07_DATA_DIR)
endif
	-rmdir $(DESTDIR)$(datarootdir)/LanlGeoMag
