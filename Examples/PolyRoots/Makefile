# Very simple makefile illustrating how to use pkg-config to compile


quartic: quartic.c
	gcc quartic.c `pkg-config --cflags --libs lgm` -o quartic

clean:
	rm quartic
