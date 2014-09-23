#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CFLAGS    = `pkg-config --cflags lgm` -Wall
CC        = gcc

FILES1    =  McIlwain_L.c
FILES2    =  McIlwain_L_2.c
FILES3    =  McIlwain_L_3.c
FILES4    =  McIlwain_L_4.c
FILES5    =  McIlwain_L_Test.c
FILES6    =  McIlwain_L_CDIP.c

all   : McIlwain_L McIlwain_L_2 McIlwain_L_3 McIlwain_L_4 McIlwain_L_Test McIlwain_L_CDIP

McIlwain_L : Makefile $(FILES1) McIlwain_L.o
	$(CC) $(FILES1) $(LIBS) $(CFLAGS) -o McIlwain_L

McIlwain_L_2 : Makefile $(FILES2) McIlwain_L_2.o
	$(CC) $(FILES2) $(LIBS) $(CFLAGS) -o McIlwain_L_2

McIlwain_L_3 : Makefile $(FILES3) McIlwain_L_3.o
	$(CC) $(FILES3) $(LIBS) $(CFLAGS) -o McIlwain_L_3

McIlwain_L_4 : Makefile $(FILES4) McIlwain_L_4.o
	$(CC) $(FILES4) $(LIBS) $(CFLAGS) -o McIlwain_L_4

McIlwain_L_Test : Makefile $(FILES5) McIlwain_L_Test.o
	$(CC) $(FILES5) $(LIBS) $(CFLAGS) -o McIlwain_L_Test

McIlwain_L_CDIP : Makefile $(FILES6) McIlwain_L_CDIP.o
	$(CC) $(FILES6) $(LIBS) $(CFLAGS) -o McIlwain_L_CDIP


clean :
	rm -f McIlwain_L McIlwain_L_2 McIlwain_L_3 McIlwain_L_4 McIlwain_L_Test McIlwain_L_CDIP *.o
