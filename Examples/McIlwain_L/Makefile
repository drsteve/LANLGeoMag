#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm`
CC        = gcc

FILES1    =  McIlwain_L.c
FILES2    =  McIlwain_L_2.c
FILES3    =  McIlwain_L_3.c

all   : McIlwain_L McIlwain_L_2 McIlwain_L_3

McIlwain_L : Makefile $(FILES1) McIlwain_L.o
	$(CC) $(FILES1) $(LIBS) $(CFLAGS) -o McIlwain_L

McIlwain_L_2 : Makefile $(FILES2) McIlwain_L_2.o
	$(CC) $(FILES2) $(LIBS) $(CFLAGS) -o McIlwain_L_2

McIlwain_L_3 : Makefile $(FILES3) McIlwain_L_3.o
	$(CC) $(FILES3) $(LIBS) $(CFLAGS) -o McIlwain_L_3


clean :
	rm McIlwain_L McIlwain_L_2 McIlwain_L_3 *.o
