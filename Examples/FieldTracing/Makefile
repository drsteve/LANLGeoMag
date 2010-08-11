# Very simple makefile illustrating how to use pkg-config to compile


SimpleTrace: SimpleTrace.c
	gcc SimpleTrace.c -mtune=core2 `pkg-config --cflags --libs lgm` -o SimpleTrace 

clean:
	rm SimpleTrace
