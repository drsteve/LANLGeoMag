#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  Elapsed.c

all   : Elapsed

Elapsed : Makefile $(FILES1) Elapsed.o
	$(CC) $(FILES1) $(LIBS) $(CFLAGS) -o Elapsed



clean :
	rm Elapsed *.o
