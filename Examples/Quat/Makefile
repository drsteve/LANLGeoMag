# Very simple makefile illustrating how to use pkg-config to compile


all: Ex1 Ex2

Ex1: Ex1.c
	gcc Ex1.c -Wall `pkg-config --cflags --libs lgm` -o Ex1


Ex2: Ex2.c
	gcc Ex2.c -Wall `pkg-config --cflags --libs lgm` -o Ex2

clean:
	rm Ex1 Ex2
