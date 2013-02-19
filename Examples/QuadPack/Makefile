# Very simple makefile illustrating how to use pkg-config to compile


slerp: ex1.c
	gcc ex1.c `pkg-config --cflags --libs lgm` -o ex1

clean:
	rm ex1
