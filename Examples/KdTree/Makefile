# Very simple makefile illustrating how to use pkg-config to compile

all: test DstNNs TLE_NNs

test: test.c
	gcc test.c `pkg-config --cflags --libs lgm` -o test

DstNNs: DstNNs.c
	gcc DstNNs.c `pkg-config --cflags --libs lgm` -o DstNNs

TLE_NNs: TLE_NNs.c
	gcc TLE_NNs.c `pkg-config --cflags --libs lgm` -o TLE_NNs

clean:
	rm test DstNNs TLE_NNs
