#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  

all   : LstarVersusPA

LstarVersusPA : Makefile $(FILES1) LstarVersusPA.o
	$(CC) $(FILES1) LstarVersusPA.o $(LIBS) $(CFLAGS) -o LstarVersusPA


clean :
	rm LstarVersusPA *.o
