# Very simple makefile illustrating how to use pkg-config to compile

HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

all: test DstNNs TLE_NNs magConj

test: test.c
	gcc test.c $(LGMFLAGS) $(HDF5FLAGS) -o test

DstNNs: DstNNs.c
	gcc DstNNs.c $(LGMFLAGS) $(HDF5FLAGS) -o DstNNs

TLE_NNs: TLE_NNs.c
	gcc TLE_NNs.c $(LGMFLAGS) $(HDF5FLAGS) -o TLE_NNs

magConj: magConj.c
	gcc magConj.c $(LGMFLAGS) $(HDF5FLAGS) -o magConj

clean:
	rm test DstNNs TLE_NNs magConj
