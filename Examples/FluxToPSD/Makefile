# Very simple makefile illustrating how to use pkg-config to compile

all: f2p f2p_v2 MuTest

MuTest: MuTest.c
	gcc MuTest.c `pkg-config --cflags --libs lgm` -o MuTest

f2p: f2p.c
	gcc f2p.c `pkg-config --cflags --libs lgm` -o f2p

f2p_v2: f2p_v2.c
	gcc f2p_v2.c DumpSvg.c `pkg-config --cflags --libs lgm` -o f2p_v2

clean:
	rm f2p
