# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: Utils.c
	gcc Utils.c `pkg-config --cflags --libs lgm` -o Utils

clean:
	rm Utils
