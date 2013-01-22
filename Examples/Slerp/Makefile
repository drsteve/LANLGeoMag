# Very simple makefile illustrating how to use pkg-config to compile


slerp: slerp.c
	gcc slerp.c `pkg-config --cflags --libs lgm` -o slerp

clean:
	rm slerp
