# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

all   : CoordQuickStart CoordTransTest compareGSE

CoordQuickStart: CoordQuickStart.c
	gcc CoordQuickStart.c -Wall $(LGMFLAGS) $(HDF5FLAGS) -o CoordQuickStart

CoordTransTest: CoordTransTest.c
	gcc CoordTransTest.c -Wall $(LGMFLAGS) $(HDF5FLAGS) -o CoordTransTest

compareGSE: compareGSE.c
	gcc compareGSE.c -Wall $(LGMFLAGS) $(HDF5FLAGS) -o compareGSE

clean:
	rm CoordQuickStart CoordTransTest compareGSE
