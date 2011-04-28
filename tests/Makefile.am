## Process this file with automake to produce Makefile.in

lgm_includes=$(top_builddir)/libLanlGeoMag/Lgm/
check_PROGRAMS = check_libLanlGeoMag check_IsoTimeStringToDateTime check_McIlwain_L
TESTS = check_libLanlGeoMag check_IsoTimeStringToDateTime check_McIlwain_L

check_libLanlGeoMag_SOURCES = check_libLanlGeoMag.c $(lgm_includes)/Lgm_CTrans.h
check_libLanlGeoMag_CFLAGS = @CHECK_CFLAGS@
check_libLanlGeoMag_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_IsoTimeStringToDateTime_SOURCES = check_IsoTimeStringToDateTime.c $(lgm_includes)/Lgm_CTrans.h
check_IsoTimeStringToDateTime_CFLAGS = @CHECK_CFLAGS@
check_IsoTimeStringToDateTime_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_McIlwain_L_SOURCES = check_McIlwain_L.c $(lgm_includes)/Lgm_CTrans.h
check_McIlwain_L_CFLAGS = @CHECK_CFLAGS@
check_McIlwain_L_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@
