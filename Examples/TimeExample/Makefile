# Very simple makefile illustrating how to use pkg-config to compile


TimeEx: TimeEx.c
	gcc TimeEx.c `pkg-config --cflags --libs lgm` -o TimeEx

clean:
	rm TimeEx
