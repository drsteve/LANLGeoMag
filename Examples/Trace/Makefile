# Very simple makefile illustrating how to use pkg-config to compile

all: MagTrace MagTrace2

MagTrace: MagTrace.c
	gcc MagTrace.c `pkg-config --cflags --libs lgm` -fopenmp -o MagTrace

MagTrace2: MagTrace2.c
	gcc MagTrace2.c `pkg-config --cflags --libs lgm` -fopenmp -o MagTrace2

clean:
	rm MagTrace MAgTrace2
