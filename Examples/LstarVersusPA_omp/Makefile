#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
HDF5FLAGS = `pkg-config hdf5 --cflags --libs`
LGMFLAGS  = `pkg-config lgm --cflags --libs`
CC        = gcc

FILES1    = LstarVersusPA.c

all   : LstarVersusPA

LstarVersusPA : $(FILES1)
	$(CC) $(FILES1) $(HDF5FLAGS) $(LGMFLAGS) -o LstarVersusPA


clean :
	rm LstarVersusPA *.o
