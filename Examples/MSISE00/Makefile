# Very simple makefile illustrating how to use pkg-config to compile

all: msis msis2

msis: msis.c
	gcc msis.c `pkg-config --cflags --libs lgm` -lgomp -o msis

msis2: msis2.c
	gcc msis2.c `pkg-config --cflags --libs lgm` -lgomp -o msis2

clean:
	rm msis msis2
