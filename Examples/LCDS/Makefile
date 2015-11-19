# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

all: LstarVersusPA LCDS

LstarVersusPA : LstarVersusPA.c
	gcc LstarVersusPA.c $(LGMFLAGS) $(HDF5FLAGS) -o LstarVersusPA

LCDS: LCDS.c
	gcc LCDS.c $(LGMFLAGS) $(HDF5FLAGS) -o LCDS

clean:
	rm LCDS LstarVersusPA
