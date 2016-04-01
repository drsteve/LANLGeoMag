# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

all: QinDenton

QinDenton: QinDenton.c
	gcc QinDenton.c $(LGMFLAGS) $(HDF5FLAGS) -o QinDenton

clean:
	rm QinDenton
