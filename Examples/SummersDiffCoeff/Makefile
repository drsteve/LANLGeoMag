# Very simple makefile illustrating how to use pkg-config to compile

all: SummersDiffCoeff SummersDiffCoeffBounceAvg

SummersDiffCoeff: SummersDiffCoeff.c
	gcc SummersDiffCoeff.c `pkg-config --cflags --libs lgm` -o SummersDiffCoeff

SummersDiffCoeffBounceAvg: SummersDiffCoeffBounceAvg.c
	gcc SummersDiffCoeffBounceAvg.c `pkg-config --cflags --libs lgm` -o SummersDiffCoeffBounceAvg


clean:
	rm SummersDiffCoeff SummersDiffCoeffBounceAvg
