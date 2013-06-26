# Very simple makefile illustrating how to use pkg-config to compile

all: test DstNNs

test: test.c
	gcc test.c `pkg-config --cflags --libs lgm` -o test

DstNNs: DstNNs.c
	gcc DstNNs.c `pkg-config --cflags --libs lgm` -o DstNNs

clean:
	rm test DstNNs
