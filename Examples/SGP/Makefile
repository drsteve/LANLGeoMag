# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs`
LGMFLAGS  = `pkg-config lgm --cflags --libs`

Geodetic: SGP.c
	gcc SGP.c $(HDF5FLAGS) $(LGMFLAGS) -Wall -o SGP

clean:
	rm SGP
