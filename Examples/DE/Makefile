# Very simple makefile illustrating how to use pkg-config to compile
HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`

all   : DEQuickStart DESunVector DE421Test

DE421Test: DE421Test.c
	gcc DE421Test.c -Wall $(LGMFLAGS) $(HDF5FLAGS) -o DE421Test

DEQuickStart: DEQuickStart.c
	gcc DEQuickStart.c -Wall $(LGMFLAGS) $(HDF5FLAGS) -o DEQuickStart

DESunVector: DESunVector.c
	gcc DESunVector.c -Wall $(LGMFLAGS) $(HDF5FLAGS) -o DESunVector

clean:
	rm DEQuickStart DESunVector DE421Test
