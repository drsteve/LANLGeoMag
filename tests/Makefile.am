## Process this file with automake to produce Makefile.in

lgm_includes=$(top_srcdir)/libLanlGeoMag/Lgm/
#check_PROGRAMS = check_libLanlGeoMag check_IsoTimeStringToDateTime check_McIlwain_L check_PolyRoots check_Sgp4
#TESTS = check_libLanlGeoMag check_IsoTimeStringToDateTime check_McIlwain_L check_PolyRoots check_Sgp4
check_PROGRAMS = check_libLanlGeoMag check_McIlwain_L check_PolyRoots check_Sgp4 check_DE421
TESTS          = check_libLanlGeoMag check_McIlwain_L check_PolyRoots check_Sgp4 check_DE421

check_libLanlGeoMag_SOURCES = check_libLanlGeoMag.c $(lgm_includes)/Lgm_CTrans.h
check_libLanlGeoMag_CFLAGS = @CHECK_CFLAGS@
check_libLanlGeoMag_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

#check_IsoTimeStringToDateTime_SOURCES = check_IsoTimeStringToDateTime.c $(lgm_includes)/Lgm_CTrans.h
#check_IsoTimeStringToDateTime_CFLAGS = @CHECK_CFLAGS@
#check_IsoTimeStringToDateTime_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_McIlwain_L_SOURCES = check_McIlwain_L.c $(lgm_includes)/Lgm_CTrans.h
check_McIlwain_L_CFLAGS = @CHECK_CFLAGS@
check_McIlwain_L_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_PolyRoots_SOURCES = check_PolyRoots.c $(lgm_includes)/Lgm_CTrans.h
check_PolyRoots_CFLAGS = @CHECK_CFLAGS@
check_PolyRoots_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_Sgp4_SOURCES = check_Sgp4.c $(lgm_includes)/Lgm_CTrans.h
check_Sgp4_CFLAGS = @CHECK_CFLAGS@
check_Sgp4_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_DE421_SOURCES = check_DE421.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_Vec.h
check_DE421_CFLAGS = @CHECK_CFLAGS@
check_DE421_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

EXTRA_DIST = check_McIlwain_L_01.expected check_McIlwain_L_02.expected check_McIlwain_L_03.expected check_McIlwain_L_04.expected check_PolyRoots_01.expected check_PolyRoots_02.expected check_PolyRoots_03.expected check_PolyRoots_04.expected check_Sgp4_01.expected testpo.421
