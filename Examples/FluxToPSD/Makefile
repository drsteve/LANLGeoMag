# Very simple makefile illustrating how to use pkg-config to compile

all: f2p MuTest

MuTest: MuTest.c
	gcc MuTest.c `pkg-config --cflags --libs lgm` -o MuTest

f2p: f2p.c
	gcc f2p.c `pkg-config --cflags --libs lgm` -fopenmp -o f2p

clean:
	rm f2p
