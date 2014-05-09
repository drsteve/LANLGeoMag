# Very simple makefile illustrating how to use pkg-config to compile


QinDenton: LCDS.c
	gcc LCDS.c `pkg-config --cflags --libs lgm` -o LCDS

clean:
	rm LCDS
