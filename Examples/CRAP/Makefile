#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  ComputeFieldQuantities2.c

all   : MagCoords

MagCoords : Makefile $(FILES1) MagCoords.o
	$(CC) $(FILES1) MagCoords.o $(LIBS) $(CFLAGS) -o MagCoords


clean :
	rm MagCoords *.o
