#BASE_DIR = /home/mgh/DREAM/Dream
#DESTDIR= /usr/local
#LGM_LIB = ./libLanlGeoMag
#LGM_INC = ./libLanlGeoMag
##OMP_INC = /usr/local/src/gcc-4.3.1-obj/x86_64-unknown-linux-gnu/libgomp/
#CFLAGS =  -I$(LGM_INC) -g -Wall -O3 -mfpmath=sse,387 -mtune=opteron -funroll-all-loops -fprefetch-loop-arrays -ffast-math -fopenmp `pkg-config --cflags gtk+-2.0` `pkg-config --cflags gtkglext-1.0`
CFLAGS =  -mfpmath=387 `pkg-config --cflags gtk+-2.0 gtkglext-1.0 lgm` 

FILES  =  support.o ReadPng.o trackball.o MakeTube.o Atmosphere.o Ellipsoid.o Sphere.o xvgifwr2.o \
		  IridiumFlare.o ComputeZPSTransMatrix.o SatSelector_Interface.o  Vds_DriftShell.c\
		  SatSelector_Callbacks.o SaveRasterFileDialog.o OpenMagEphemFile.o

LIBS   = -g -lGLEW `pkg-config --cflags --libs gtk+-2.0 gthread-2.0 gtkglext-1.0 lgm`
CC     = gcc 

all: ViewDriftShell

ViewDriftShell : Makefile ViewDriftShell.o $(FILES)
	$(CC) ViewDriftShell.o $(FILES) $(LIBS) $(CFLAGS) -o ViewDriftShell `pkg-config --cflags --libs gtk+-2.0 lgm` 



clean : 
	rm *.o 

install:  
