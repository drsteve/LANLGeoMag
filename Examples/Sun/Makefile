# Very simple makefile illustrating how to use pkg-config to compile


foo: foo.c
	gcc foo.c `pkg-config --cflags --libs lgm` -o foo

clean:
	rm foo
