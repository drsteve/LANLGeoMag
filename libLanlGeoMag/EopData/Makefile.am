AM_CPPFLAGS = -DDATADIR='"$(datadir)"'

EOPDATAdir   = $(datadir)/LanlGeoMag/EopData
EOPDATA_DATA = GetEopFiles LgmEop.dat eopc04.62-now finals.daily \
               finals2000A.all finals.all finals2000A.daily tai-utc.dat \
               eopc04_IAU2000.62-now GetEopFiles.in Lgm_LeapSecondDates.dat

EXTRA_DIST = $(EOPDATA_DATA)

#
# Data install hooks are guaranteed to be run after the data is actually
# installed (Dont use install-exec-hook -- that one isnt guaranteed to run
# after the data install.)
#	if [ $(ENABLE_EOP_CRON) -eq 1 ]; then \
#
ENABLE_EOP_CRON = @ENABLE_EOP_CRON@
DAILY_CRON_DIR  = @DAILY_CRON_DIR@
install-data-hook:
	chmod +x $(DESTDIR)$(EOPDATAdir)/GetEopFiles;
#	$(EOPDATAdir)/GetEopFiles 
	if [ $(ENABLE_EOP_CRON) -eq 1 ]; then \
	ln -sf $(DESTDIR)$(EOPDATAdir)/GetEopFiles $(DAILY_CRON_DIR)/GetEopFiles;\
	fi;

uninstall-hook:
	-rmdir $(DESTDIR)$(EOPDATAdir)
