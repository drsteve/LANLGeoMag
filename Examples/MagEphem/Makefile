#CFLAGS       = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS     =  -fno-automatic 
SPICE_LIBS = '/home/mgh/SPICE/cspice/lib/cspice.a'
SPICE_INC  = '-I/home/mgh/SPICE/cspice/include'
#HDF5_LIBS  = '-L/packages/lib/hdf5-1.8.7'
#HDF5_INC   = '-I/packages/include/hdf5-1.8'
LIBS       = `pkg-config --libs lgm`
CFLAGS     = `pkg-config --cflags lgm` $(SPICE_INC)
CC           = gcc


all   : MagEphemFromFile MagEphemFromSpiceKernel

MagEphemFromFile : Makefile  MagEphemFromFile.o
	$(CC) $(FILES1) MagEphemFromFile.o $(HDF5_LIBS) $(LIBS) $(CFLAGS) -o MagEphemFromFile

MagEphemFromSpiceKernel : Makefile  MagEphemFromSpiceKernel.o
	$(CC) $(FILES1) MagEphemFromSpiceKernel.o $(HDF5_LIBS) $(LIBS) $(CFLAGS) $(SPICE_LIBS) -o MagEphemFromSpiceKernel



clean :
	rm MagEphemFromFile MagEphemFromSpiceKernel *.o
