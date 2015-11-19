# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

AE8_AP8: AE8_AP8.c
	gcc AE8_AP8.c $(LGMFLAGS) $(HDF5FLAGS) -o AE8_AP8

clean:
	rm AE8_AP8
