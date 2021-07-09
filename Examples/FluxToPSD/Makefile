# Very simple makefile illustrating how to use pkg-config to compile

all: f2p f2p_v2 MuTest FluxToPSD PSD_Versus_EA

MuTest: MuTest.c
	gcc MuTest.c `pkg-config --cflags --libs lgm` -o MuTest

f2p: f2p.c
	gcc f2p.c `pkg-config --cflags --libs lgm` -o f2p

f2p_v2: f2p_v2.c DumpSvg.c
	gcc f2p_v2.c DumpSvg.c `pkg-config --cflags --libs lgm` -o f2p_v2

FluxToPSD: FluxToPSD.c DumpSvg.c
	h5cc FluxToPSD.c DumpSvg.c CDFUtils.c `pkg-config --cflags --libs lgm` -lcdf -o FluxToPSD

clean:
	-rm -f f2p f2p_v2 FluxToPSD
