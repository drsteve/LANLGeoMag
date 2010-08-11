# Very simple makefile illustrating how to use pkg-config to compile


MultiCoord: MultiCoord.c
	gcc MultiCoord.c `pkg-config --cflags --libs lgm` -o MultiCoord

clean:
	rm MultiCoord
