# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

SimpleTrace: SimpleTrace.c
	gcc SimpleTrace.c -mtune=core2 $(LGMFLAGS) $(HDF5FLAGS) -o SimpleTrace 

AlphaOfK: AlphaOfK.c
	gcc AlphaOfK.c -mtune=core2 $(LGMFLAGS) $(HDF5FLAGS) -o AlphaOfK 

clean:
	rm SimpleTrace AlphaOfK
