# Very simple makefile illustrating how to use pkg-config to compile


CoordQuickStart: CoordQuickStart.c
	gcc CoordQuickStart.c -Wall `pkg-config --cflags --libs lgm` -o CoordQuickStart

clean:
	rm CoordQuickStart
