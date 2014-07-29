# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: Utils.c Makefile
	gcc Utils.c `pkg-config --cflags --libs lgm` -Wall -o Utils

clean:
	rm Utils
