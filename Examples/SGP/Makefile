# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: SGP.c
	gcc SGP.c `pkg-config --cflags --libs lgm` -Wall -o SGP

clean:
	rm SGP
