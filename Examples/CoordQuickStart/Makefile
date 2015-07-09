# Very simple makefile illustrating how to use pkg-config to compile

all   : CoordQuickStart CoordTransTest

CoordQuickStart: CoordQuickStart.c
	gcc CoordQuickStart.c -Wall `pkg-config --cflags --libs lgm` -o CoordQuickStart

CoordTransTest: CoordTransTest.c
	gcc CoordTransTest.c -Wall `pkg-config --cflags --libs lgm` -o CoordTransTest

clean:
	rm CoordQuickStart CoordTransTest
