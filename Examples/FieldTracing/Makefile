# Very simple makefile illustrating how to use pkg-config to compile


SimpleTrace: SimpleTrace.c
	gcc SimpleTrace.c -mtune=core2 `pkg-config --cflags --libs lgm` -o SimpleTrace 

AlphaOfK: AlphaOfK.c
	gcc AlphaOfK.c -mtune=core2 `pkg-config --cflags --libs lgm` -o AlphaOfK 

clean:
	rm SimpleTrace AlphaOfK
