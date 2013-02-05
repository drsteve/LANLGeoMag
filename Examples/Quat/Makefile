# Very simple makefile illustrating how to use pkg-config to compile


slerp: Ex1.c
	gcc Ex1.c -Wall `pkg-config --cflags --libs lgm` -o Ex1

clean:
	rm Ex1
