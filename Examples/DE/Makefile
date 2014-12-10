# Very simple makefile illustrating how to use pkg-config to compile


DEQuickStart: DEQuickStart.c
	gcc DEQuickStart.c -Wall `pkg-config --cflags --libs lgm` -o DEQuickStart

clean:
	rm DEQuickStart
