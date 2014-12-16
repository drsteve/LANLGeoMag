# Very simple makefile illustrating how to use pkg-config to compile

all   : DEQuickStart DESunVector DE421Test

DE421Test: DE421Test.c
	gcc DE421Test.c -Wall `pkg-config --cflags --libs lgm` -o DE421Test

DEQuickStart: DEQuickStart.c
	gcc DEQuickStart.c -Wall `pkg-config --cflags --libs lgm` -o DEQuickStart

DESunVector: DESunVector.c
	gcc DESunVector.c -Wall `pkg-config --cflags --libs lgm` -o DESunVector

clean:
	rm DEQuickStart DESunVector DE421Test
