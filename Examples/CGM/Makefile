# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

all: CgmCoords

CgmCoords: CgmCoords.c
	gcc CgmCoords.c $(LGMFLAGS) $(HDF5FLAGS) -o CgmCoords


clean:
	rm CgmCoords
