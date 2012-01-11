AM_CPPFLAGS = -DDATADIR='"$(datadir)"'

EOPDATAdir   = @datadir@/LanlGeoMag/EopData
EOPDATA_DATA = GetEopFiles LgmEop.dat eopc04.62-now finals.daily \
               finals2000A.all finals.all finals2000A.daily tai-utc.dat \
               eopc04_IAU2000.62-now GetEopFiles.in Lgm_LeapSecondDates.dat

#
# Data install hooks are guaranteed to be run after the data is actually
# installed (Dont use install-exec-hook -- that one isnt guaranteed to run
# after the data install.)
#
install-data-hook:
	chmod +x $(EOPDATAdir)/GetEopFiles;
	$(EOPDATAdir)/GetEopFiles 
