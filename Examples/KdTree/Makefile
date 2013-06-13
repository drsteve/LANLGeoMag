# Very simple makefile illustrating how to use pkg-config to compile


test: test.c
	gcc test.c `pkg-config --cflags --libs lgm` -o test

clean:
	rm test
