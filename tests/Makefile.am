## Process this file with automake to produce Makefile.in

lgm_includes=$(top_srcdir)/libLanlGeoMag/Lgm/
check_PROGRAMS = check_libLanlGeoMag check_ClosedField check_McIlwain_L check_PolyRoots check_Magmodels check_Sgp4 check_DE421 check_CoordTrans check_IsoTimeStringToDateTime check_Lstar check_string_funcs
TESTS          = check_libLanlGeoMag check_ClosedField check_McIlwain_L check_PolyRoots check_Magmodels check_Sgp4 check_DE421 check_CoordTrans check_IsoTimeStringToDateTime check_Lstar check_string_funcs

check_libLanlGeoMag_SOURCES = check_libLanlGeoMag.c $(lgm_includes)/Lgm_CTrans.h
check_libLanlGeoMag_CFLAGS = @CHECK_CFLAGS@
check_libLanlGeoMag_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_ClosedField_SOURCES = check_ClosedField.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_MagModelInfo.h
check_ClosedField_CFLAGS = @CHECK_CFLAGS@
check_ClosedField_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @PERL_LDFLAGS@ @CHECK_LIBS@

check_IsoTimeStringToDateTime_SOURCES = check_IsoTimeStringToDateTime.c $(lgm_includes)/Lgm_CTrans.h
check_IsoTimeStringToDateTime_CFLAGS = @CHECK_CFLAGS@
check_IsoTimeStringToDateTime_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @PERL_LDFLAGS@ @CHECK_LIBS@

check_McIlwain_L_SOURCES = check_McIlwain_L.c $(lgm_includes)/Lgm_CTrans.h
check_McIlwain_L_CFLAGS = @CHECK_CFLAGS@
check_McIlwain_L_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_Lstar_SOURCES = check_Lstar.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_LstarInfo.h $(lgm_includes)/Lgm_MagEphemInfo.h
check_Lstar_CFLAGS = @CHECK_CFLAGS@ -fopenmp
check_Lstar_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @PERL_LDFLAGS@ @CHECK_LIBS@

check_PolyRoots_SOURCES = check_PolyRoots.c $(lgm_includes)/Lgm_CTrans.h
check_PolyRoots_CFLAGS = @CHECK_CFLAGS@
check_PolyRoots_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_Sgp4_SOURCES = check_Sgp4.c $(lgm_includes)/Lgm_CTrans.h
check_Sgp4_CFLAGS = @CHECK_CFLAGS@
check_Sgp4_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_Magmodels_SOURCES = check_Magmodels.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_MagModelInfo.h
check_Magmodels_CFLAGS = @CHECK_CFLAGS@
check_Magmodels_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @PERL_LDFLAGS@ @CHECK_LIBS@

check_DE421_SOURCES = check_DE421.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_Vec.h
check_DE421_CFLAGS = @CHECK_CFLAGS@
check_DE421_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @CHECK_LIBS@

check_CoordTrans_SOURCES = check_CoordTrans.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_Vec.h
check_CoordTrans_CFLAGS = @CHECK_CFLAGS@
check_CoordTrans_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @PERL_LDFLAGS@ @CHECK_LIBS@

check_string_funcs_SOURCES = check_string_funcs.c $(lgm_includes)/Lgm_CTrans.h $(lgm_includes)/Lgm_Misc.h
check_string_funcs_CFLAGS = @CHECK_CFLAGS@
check_string_funcs_LDADD = $(top_builddir)/libLanlGeoMag/.libs/libLanlGeoMag.a @PERL_LDFLAGS@ @CHECK_LIBS@

EXTRA_DIST = check_Lstar.expected check_McIlwain_L_01.expected check_McIlwain_L_02.expected check_McIlwain_L_03.expected check_McIlwain_L_04.expected check_McIlwain_L_05.expected check_McIlwain_L_06.expected check_McIlwain_L_07.expected check_McIlwain_L_08.expected check_PolyRoots_01.expected check_PolyRoots_02.expected check_PolyRoots_03.expected check_PolyRoots_04.expected check_Sgp4_01.expected testpo.421 check_CoordTrans.expected check_CoordTransNoEph.expected check_CoordDipoleTilt.expected check_Magmodels_01.expected check_ClosedField_01.expected
