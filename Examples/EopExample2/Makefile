# Very simple makefile illustrating how to use pkg-config to compile


EopParams: EopParams.c
	gcc EopParams.c `pkg-config --cflags --libs lgm` -o EopParams

clean:
	rm EopParams
