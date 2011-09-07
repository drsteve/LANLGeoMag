#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  

all   : MagEphemFromTLE MagEphemFromFile MagEphemFromTLE2 MagEphemFromPosition MagEphemFromPosition2

MagEphemFromFile : Makefile $(FILES1) MagEphemFromFile.o
	$(CC) $(FILES1) MagEphemFromFile.o $(LIBS) $(CFLAGS) -o MagEphemFromFile

MagEphemFromTLE : Makefile $(FILES1) MagEphemFromTLE.o
	$(CC) $(FILES1) MagEphemFromTLE.o $(LIBS) $(CFLAGS) -o MagEphemFromTLE

MagEphemFromTLE2 : Makefile $(FILES1) MagEphemFromTLE2.o
	$(CC) $(FILES1) MagEphemFromTLE2.o $(LIBS) $(CFLAGS) -o MagEphemFromTLE2

MagEphemFromPosition : Makefile $(FILES1) MagEphemFromPosition.o
	$(CC) $(FILES1) MagEphemFromPosition.o $(LIBS) $(CFLAGS) -o MagEphemFromPosition

MagEphemFromPosition2 : Makefile $(FILES1) MagEphemFromPosition2.o
	$(CC) $(FILES1) MagEphemFromPosition2.o $(LIBS) $(CFLAGS) -o MagEphemFromPosition2


clean :
	rm MagEphemFromTLE MagEphemFromFile MagEphemFromTLE2 MagEphemFromPosition MagEphemFromPosition2 *.o
