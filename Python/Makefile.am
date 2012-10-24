#Don't forget: probably need to do stuff with EXTRA_DIST
LGMSRCDIR=../libLanlGeoMag/

all-local: $(LGMSRCDIR)/libLanlGeoMag.la
	LIBDIR=$(libdir) $(PYTHON) setup.py build

install-exec-local:
	if [ -n "$(DESTDIR)" ]; then LIBDIR=$(libdir) $(PYTHON) setup.py install --root=$(DESTDIR) --prefix=$(prefix); else LIBDIR=$(libdir) $(PYTHON) setup.py install --prefix=$(prefix); fi

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/lgmpy*

#http://www.gnu.org/software/automake/manual/html_node/Clean.html
mostlyclean-local: mostlyclean-local-check

.PHONY: mostlyclean-local-check

mostlyclean-local-check:
	-rm -rf build
