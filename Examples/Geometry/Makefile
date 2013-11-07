# Very simple makefile illustrating how to use pkg-config to compile


Geodetic: SphereRay.c
	gcc SphereRay.c -Wall -O2 `pkg-config --cflags --libs lgm` -o SphereRay

clean:
	rm SphereRay
