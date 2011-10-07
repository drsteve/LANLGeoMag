SUBDIRS = Lgm

AM_CPPFLAGS = -DDATADIR='"$(datadir)"'

# Location to keep Eop-related stuff
LGM_EOP_DATA_DIR   = @prefix@/share/LanlGeoMag/EopData
LGM_INDEX_DATA_DIR = @prefix@/share/LanlGeoMag/Data
includedir         = @prefix@/include/Lgm
#libdir             = @prefix@/lib
pkgconfigdir       = $(libdir)/pkgconfig
pkgconfig_DATA     = lgm.pc

if LINUX
    CRONDIR = /etc/cron.daily
else
    CRONDIR = /etc/periodic/daily
endif

AM_CFLAGS = -fopenmp -lgomp

lib_LTLIBRARIES          = libLanlGeoMag.la
libLanlGeoMag_la_SOURCES =  Lgm_AlphaOfK.c B_FromScatteredData.c ComputeLstar.c DriftShell.c IntegralInvariant.c LFromIBmM.c \
	                        Lgm_B_internal.c Lgm_CTrans.c Lgm_DateAndTime.c Lgm_Eop.c Lgm_IGRF.c Lgm_InitMagInfo.c \
                            Lgm_MaxwellJuttner.c Lgm_Nutation.c Lgm_Octree.c Lgm_Quat.c Lgm_Sgp.c Lgm_SimplifiedMead.c  Lgm_SunPosition.c \
                            Lgm_Trace.c Lgm_TraceToEarth.c Lgm_TraceToSphericalEarth.c Lgm_Vec.c MagStep.c Lgm_QuadPack3.c \
                            Lgm_QuadPack.c Lgm_Cgm.c quicksort.c SbIntegral.c T87.c T89.c TraceLine.c Lgm_TraceToMinBSurf.c  \
                            TraceToMinRdotB.c Lgm_TraceToMirrorPoint.c TraceToSMEquat.c T01S.c Tsyg_T01s.c TS04.c Tsyg2004.c \
                            TS04_opt.c Tsyg2004_Optimized.c Lgm_PolyRoots.c Lgm_SummersDiffCoeff.c \
                            W.c Lgm_InitMagEphemInfo.c Lgm_AE8_AP8.c OP77.c OlsenPfitzerStatic.c IsoTimeStringToDateTime.c \
			                size.c Lgm_FluxToPsd.c xvgifwr2.c praxis.c Lgm_SphHarm.c Lgm_McIlwain_L.c Lgm_ElapsedTime.c \
			                Lgm_ComputeLstarVersusPA.c Lgm_MagEphemWrite.c brent.c Lgm_CdipMirrorLat.c ComputeI_FromMltMlat.c \
                            Lgm_QinDenton.c Lgm_DiffCoeff_param.c Lgm_AE_index.c Lgm_Misc.c Lgm_HDF5.c Lgm_GradB.c Lgm_VelStep.c #Lgm_Enkf_1D_RadDiff.c #Lgm_FastPowPoly.c
                            
include_HEADERS = Lgm/*.h

EOP_DATA_FILES = EopData/GetEopFiles EopData/LgmEop.dat EopData/eopc04.62-now \
	EopData/finals.daily EopData/finals2000A.all EopData/finals.all \
	EopData/finals2000A.daily EopData/tai-utc.dat \
	EopData/eopc04_IAU2000.62-now EopData/GetEopFiles.in \
	EopData/Lgm_LeapSecondDates.dat

QINDENTON_DIR = QinDenton

#EXTRA_DIST = EopData
EXTRA_DIST = $(EOP_DATA_FILES) eop_data_file_list $(QINDENTON_DIR)





install-exec-hook:
	mkdir -p $(LGM_EOP_DATA_DIR)
#	sed -e "s+CHANGE_ME+$(LGM_EOP_DATA_DIR)+" EopData/GetEopFiles.in > EopData/GetEopFiles
	chmod +x EopData/GetEopFiles
	cp EopData/GetEopFiles $(LGM_EOP_DATA_DIR)/.
	cp EopData/LgmEop.dat $(LGM_EOP_DATA_DIR)/.
	cp EopData/Lgm_LeapSecondDates.dat $(LGM_EOP_DATA_DIR)/.
	cp EopData/eopc04.62-now $(LGM_EOP_DATA_DIR)/.
	cp EopData/eopc04_IAU2000.62-now $(LGM_EOP_DATA_DIR)/.
	cp EopData/finals2000A.all $(LGM_EOP_DATA_DIR)/.
	cp EopData/finals2000A.daily $(LGM_EOP_DATA_DIR)/.
	cp EopData/finals.all $(LGM_EOP_DATA_DIR)/.
	cp EopData/finals.daily $(LGM_EOP_DATA_DIR)/.
	cp EopData/tai-utc.dat $(LGM_EOP_DATA_DIR)/.
	mkdir -p $(LGM_INDEX_DATA_DIR)
	rsync -ar QinDenton $(LGM_INDEX_DATA_DIR)/.
#	test -w "$(CRONDIR)" && ln -sf $(LGM_EOP_DATA_DIR)/GetEopFiles $(CRONDIR)/GetEopFiles
