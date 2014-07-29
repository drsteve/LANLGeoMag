# Very simple makefile illustrating how to use pkg-config to compile


slerp: slerp.c Makefile
	gcc slerp.c `pkg-config --cflags --libs lgm` -Wall -o slerp

clean:
	rm slerp
