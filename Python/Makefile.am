#Don't forget: probably need to do stuff with EXTRA_DIST
LGMSRCDIR=../libLanlGeoMag/

all-local: $(LGMSRCDIR)/libLanlGeoMag.la
	$(PYTHON) setup.py build

install-exec-local:
	if [ -n $(DESTDIR) ]; then
		$(PYTHON) setup.py install --root=$(DESTDIR) --prefix=$(prefix)
	else
		$(PYTHON) setup.py install --prefix=$(prefix)
	fi

uninstall-local:
	rm -f $(DESTDIR)$(libdir)/python*/site-packages/lgmpy*

#http://www.gnu.org/software/automake/manual/html_node/Clean.html
mostlyclean-local: mostlyclean-local-check

.PHONY: mostlyclean-local-check

mostlyclean-local-check:
	-rm -rf build
