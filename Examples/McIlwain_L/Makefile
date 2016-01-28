#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
HDF5FLAGS = `pkg-config hdf5 --libs --cflags`
LGMFLAGS  = `pkg-config lgm --libs --cflags` -Wall
CC        = gcc

FILES1    =  McIlwain_L.c
FILES2    =  McIlwain_L_2.c
FILES3    =  McIlwain_L_3.c
FILES4    =  McIlwain_L_4.c
FILES5    =  McIlwain_L_Test.c

all   : McIlwain_L McIlwain_L_2 McIlwain_L_3 McIlwain_L_4 McIlwain_L_Test

McIlwain_L : $(FILES1)
	$(CC) $(FILES1) $(HDF5FLAGS) $(LGMFLAGS) -o McIlwain_L

McIlwain_L_2 : $(FILES2)
	$(CC) $(FILES2) $(HDF5FLAGS) $(LGMFLAGS) -o McIlwain_L_2

McIlwain_L_3 : $(FILES3)
	$(CC) $(FILES3) $(HDF5FLAGS) $(LGMFLAGS) -o McIlwain_L_3

McIlwain_L_4 : $(FILES4)
	$(CC) $(FILES4) $(HDF5FLAGS) $(LGMFLAGS) -o McIlwain_L_4

McIlwain_L_Test : $(FILES5)
	$(CC) $(FILES5) $(HDF5FLAGS) $(LGMFLAGS) -o McIlwain_L_Test


clean :
	rm -f McIlwain_L McIlwain_L_2 McIlwain_L_3 McIlwain_L_4 McIlwain_L_Test *.o
