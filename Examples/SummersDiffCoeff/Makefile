# Very simple makefile illustrating how to use pkg-config to compile


SummersDiffCoeff: SummersDiffCoeff.c
	gcc SummersDiffCoeff.c `pkg-config --cflags --libs lgm` -o SummersDiffCoeff

clean:
	rm SummersDiffCoeff
