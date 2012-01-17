#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  

all   : Lstar 

Lstar : Makefile $(FILES1) Lstar.o
	$(CC) $(FILES1) Lstar.o $(LIBS) $(CFLAGS) -o Lstar



clean :
	rm Lstar *.o
