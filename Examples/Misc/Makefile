# Very simple makefile illustrating how to use pkg-config to compile


All: Misc.c
	gcc Misc.c `pkg-config --cflags --libs lgm` -o Misc

clean:
	rm Misc
