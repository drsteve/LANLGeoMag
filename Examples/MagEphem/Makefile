#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  ComputeLstarVersusPA.c MagEphemWrite.o

all   : MagEphemFromTLE MagEphemFromFile MagEphemFromTLE2

MagEphemFromFile : Makefile $(FILES1) MagEphemFromFile.o
	$(CC) $(FILES1) MagEphemFromFile.o $(LIBS) $(CFLAGS) -o MagEphemFromFile

MagEphemFromTLE : Makefile $(FILES1) MagEphemFromTLE.o
	$(CC) $(FILES1) MagEphemFromTLE.o $(LIBS) $(CFLAGS) -o MagEphemFromTLE

MagEphemFromTLE2 : Makefile $(FILES1) MagEphemFromTLE2.o
	$(CC) $(FILES1) MagEphemFromTLE2.o $(LIBS) $(CFLAGS) -o MagEphemFromTLE2


clean :
	rm MagEphemFromTLE *.o
