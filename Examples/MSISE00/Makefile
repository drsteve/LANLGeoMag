# Very simple makefile illustrating how to use pkg-config to compile


msis: msis.c
	gcc msis.c `pkg-config --cflags --libs lgm` -o msis

clean:
	rm msis
