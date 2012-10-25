#Don't forget: probably need to do stuff with EXTRA_DIST
LGMSRCDIR=../libLanlGeoMag/
lgmpydir = $(pythondir)/lgmpy
lgmpy_PYTHON = lgmpy/*.py
nodist_lgmpy_PYTHON = Lgm_Wrap.py

all-local: Lgm_Wrap.py

Lgm_Wrap.py: $(LGMSRCDIR)/libLanlGeoMag.la $(LGMSRCDIR)/Lgm/*.h Makefile
	ctypesgen.py -lLanlGeoMag --no-macro-warnings -o Lgm_Wrap.py --compile-libdir=$(LGMSRCDIR)/.libs --runtime-libdir=$(libdir) --includedir=$(LGMSRCDIR)/Lgm --includedir=$(LGMSRCDIR) $(LGMSRCDIR)/Lgm/*.h >ctypesgen.log 2>ctypesgen.log

uninstall-local:
	-rm -rf $(lgmpydir)

#http://www.gnu.org/software/automake/manual/html_node/Clean.html
mostlyclean-local: mostlyclean-local-check

.PHONY: mostlyclean-local-check

mostlyclean-local-check:
	-rm -f Lgm_Wrap.py ctypesgen.log
