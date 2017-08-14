# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: EllipsoidRay.c
	gcc EllipsoidRay.c -Wall -O2 `pkg-config --cflags --libs lgm` -o EllipsoidRay

clean:
	rm EllipsoidRay
