SUBDIRS = Lgm EopData

# Directory where Eop-related stuff will be installed
LGM_EOP_DATA_DIR   = ${prefix}/share/LanlGeoMag/EopData

# Directory where TS07-related files will be installed
LGM_TS07_DATA_DIR = ${prefix}/share/LanlGeoMag/Data

# Directory where Index-related stuff will be installed
LGM_INDEX_DATA_DIR = ${prefix}/share/LanlGeoMag/Data

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
                            Lgm_PolyRoots.c Lgm_SummersDiffCoeff.c Lgm_B_Dungey.c Lgm_FastPowPoly.c Tsyg2007.c TS07.c\
                            W.c Lgm_InitMagEphemInfo.c Lgm_AE8_AP8.c OP77.c OlsenPfitzerStatic.c IsoTimeStringToDateTime.c \
			                size.c Lgm_FluxToPsd.c xvgifwr2.c praxis.c Lgm_SphHarm.c Lgm_McIlwain_L.c Lgm_ElapsedTime.c \
			                Lgm_ComputeLstarVersusPA.c Lgm_MagEphemWrite.c brent.c Lgm_CdipMirrorLat.c ComputeI_FromMltMlat.c \
                            Lgm_QinDenton.c Lgm_DiffCoeff_param.c Lgm_AE_index.c Lgm_Misc.c Lgm_HDF5.c Lgm_GradB.c Lgm_VelStep.c Lgm_Utils.c #Lgm_Enkf_1D_RadDiff.c #Lgm_FastPowPoly.c
libLanlGeoMag_la_LIBADD = @PERL_LDFLAGS@ -lgomp
libLanlGeoMag_la_CFLAGS = $(AM_CFLAGS) @PERL_CFLAGS@ -fopenmp
libLanlGeoMag_la_CPPFLAGS = $(AM_CPPFLAGS) -DDATADIR=\"$(datadir)\" -DHASH_FUNCTION=HASH_SFH -DHASH_BLOOM=32
                            


EXTRA_DIST = $(EOP_DATA_FILES) QinDenton








install-exec-hook:
#	chmod +x $(EOPDATAdir)/GetEopFiles; $(EOPDATAdir)/GetEopFiles
	mkdir -p $(DESTDIR)$(LGM_INDEX_DATA_DIR)
	rsync -ar QinDenton $(DESTDIR)$(LGM_INDEX_DATA_DIR)/.
	rsync -ar TS07D_FILES $(DESTDIR)$(LGM_TS07_DATA_DIR)/.
