CFLAGS    = `pkg-config --cflags lgm`
#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --libs lgm`
CC        = gcc

FILES1    =  

all   : TestValladoEx315

TestValladoEx315 : Makefile $(FILES1) TestValladoEx315.o
	$(CC) $(FILES1) TestValladoEx315.o $(LIBS) -o TestValladoEx315

clean :
	rm TestValladoEx315 *.o
