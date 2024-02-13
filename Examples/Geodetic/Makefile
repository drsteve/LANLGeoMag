# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: Geodetic.c
	gcc Geodetic.c `pkg-config --cflags --libs lgm` -o Geodetic

Geodetic2: Geodetic2.c
	gcc Geodetic2.c `pkg-config --cflags --libs lgm` -o Geodetic2

clean:
	rm Geodetic Geodetic2
