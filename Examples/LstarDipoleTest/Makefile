#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LGM      = `pkg-config lgm --cflags --libs`
HDF5      = `pkg-config hdf5 --cflags --libs &2> /dev/null`
CC        = gcc

all   : LstarVersusPA

LstarVersusPA : LstarVersusPA.c
	$(CC) ComputeLstarVersusPA.c LstarVersusPA.c $(HDF5) $(LGM) -o LstarVersusPA


clean :
	rm LstarVersusPA *.o
