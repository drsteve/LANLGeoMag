AM_CPPFLAGS = -DDATADIR='"$(datadir)"'

pkgincludedir      = @includedir@/Lgm
pkginclude_HEADERS =        Lgm_CTrans.h Lgm_Eop.h Lgm_FieldIntInfo.h Lgm_IGRF.h Lgm_LstarInfo.h \
                            Lgm_MagModelInfo.h Lgm_Octree.h Lgm_QuadPack.h Lgm_Quat.h Lgm_Sgp.h Lgm_Vec.h Lgm_WGS84.h  \
                            Lgm_MagEphemInfo.h Lgm_AE8_AP8.h Lgm_DynamicMemory.h Lgm_FluxToPsd.h size.h Lgm_MaxwellJuttner.h \
                            Lgm_SphHarm.h quicksort.h  Lgm_ElapsedTime.h Lgm_PolyRoots.h Lgm_SummersDiffCoeff.h \
                            Lgm_QinDenton.h Lgm_FastPowPoly.h Lgm_Misc.h Lgm_Constants.h Lgm_DFI_RBF.h uthash.h \
                            Lgm_HDF5.h Lgm_AE_index.h qsort.h Lgm_Tsyg2004.h Lgm_Utils.h Lgm_Tsyg2007.h Lgm_Metadata.h \
                            Lgm_Tsyg1996.h Lgm_Tsyg2001.h


uninstall-hook:
	-rmdir $(pkgincludedir)
