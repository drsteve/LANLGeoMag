#Don't forget: probably need to do stuff with EXTRA_DIST
LGMSRCDIR = $(top_srcdir)/libLanlGeoMag/
lgmpydir = $(pythondir)/lgmpy
lgmpy_PYTHON = $(srcdir)/lgmpy/*.py
nodist_lgmpy_PYTHON = Lgm_Wrap.py
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`

all-local: Lgm_Wrap_gen.py #make generic wrapper at build time

Lgm_Wrap_gen.py: $(top_builddir)/libLanlGeoMag/libLanlGeoMag.la $(LGMSRCDIR)/Lgm/*.h
	PYTHONPATH=$(srcdir)/ctypesgen:${PYTHONPATH} $(PYTHON) $(srcdir)/ctypesgen/ctypesgen.py $(HDF5FLAGS) -lLanlGeoMag --no-macro-warnings -o Lgm_Wrap_gen.py --compile-libdir=$(top_builddir)/libLanlGeoMag/.libs --runtime-libdir=RUNTIME_LIBDIR --includedir=$(LGMSRCDIR)/Lgm --includedir=$(LGMSRCDIR) $(LGMSRCDIR)/Lgm/*.h >ctypesgen.log 2>ctypesgen.log
	-rm -f $(srcdir)/ctypesgen/ctypesgencore/*.pyc
	-rm -f $(srcdir)/ctypesgen/ctypesgencore/*/*.pyc

#only make final Lgm_Wrap on install, so it can rewrite the library directory
Lgm_Wrap.py: FORCE Lgm_Wrap_gen.py
	$(SED) -e 's|add_library_search_dirs(\['\''RUNTIME_LIBDIR'\''\])|add_library_search_dirs(\['\''$(libdir)'\''\])|' Lgm_Wrap_gen.py > Lgm_Wrap.py

uninstall-hook:
	-rmdir $(lgmpydir)

#http://www.gnu.org/software/automake/manual/html_node/Clean.html
mostlyclean-local: mostlyclean-local-check

FORCE:

.PHONY: mostlyclean-local-check FORCE

mostlyclean-local-check:
	-rm -f Lgm_Wrap.py Lgm_Wrap_gen.py ctypesgen.log
