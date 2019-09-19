# Very simple makefile illustrating how to use pkg-config to compile

all: SimpleTrace 

SimpleTrace: SimpleTrace.c
	gcc SimpleTrace.c `pkg-config --cflags --libs lgm` -o SimpleTrace


clean:
	rm SimpleTrace 
