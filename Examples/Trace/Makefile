# Very simple makefile illustrating how to use pkg-config to compile


MagTrace: MagTrace.c
	gcc MagTrace.c `pkg-config --cflags --libs lgm` -fopenmp -o MagTrace

clean:
	rm MagTrace
