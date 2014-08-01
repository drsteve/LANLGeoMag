# Very simple makefile illustrating how to use pkg-config to compile


quartic: quartic.c Makefile
	gcc quartic.c `pkg-config --cflags --libs lgm` -Wall -o quartic

clean:
	rm quartic
