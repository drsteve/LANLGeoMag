#CFLAGS       = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS     =  -fno-automatic 
SPICE_LIBS = '/home/mgh/SPICE/cspice/lib/cspice.a'
SPICE_INC  = '-I/home/mgh/SPICE/cspice/include'
LIBS       = `pkg-config --libs lgm`
CFLAGS     = `pkg-config --cflags lgm` $(SPICE_INC)
CC           = gcc


all   : MagEphemFromFile MagEphemFromSpiceKernel

MagEphemFromFile : Makefile  MagEphemFromFile.o
	$(CC) $(FILES1) MagEphemFromFile.o $(LIBS) $(CFLAGS) -o MagEphemFromFile

MagEphemFromSpiceKernel : Makefile  MagEphemFromSpiceKernel.o
	$(CC) $(FILES1) MagEphemFromSpiceKernel.o $(LIBS) $(CFLAGS) $(SPICE_LIBS) -o MagEphemFromSpiceKernel



clean :
	rm MagEphemFromFile MagEphemFromSpiceKernel *.o
