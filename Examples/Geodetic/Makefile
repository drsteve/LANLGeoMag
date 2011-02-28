# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: Geodetic.c
	gcc Geodetic.c `pkg-config --cflags --libs lgm` -o Geodetic

clean:
	rm Geodetic
