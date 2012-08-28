# Very simple makefile illustrating how to use pkg-config to compile

all: IsoTime IsoTime2

IsoTime: IsoTime.c
	gcc IsoTime.c `pkg-config --cflags --libs lgm` -o IsoTime

IsoTime2: IsoTime2.c
	gcc IsoTime2.c `pkg-config --cflags --libs lgm` -o IsoTime2

clean:
	rm IsoTime IsoTime2
