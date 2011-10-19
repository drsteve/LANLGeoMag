# Very simple makefile illustrating how to use pkg-config to compile

all: SummersDiffCoeff SummersDiffCoeffBounceAvg SummersDiffCoeffDerivBounceAvg

SummersDiffCoeff: SummersDiffCoeff.c
	gcc SummersDiffCoeff.c `pkg-config --cflags --libs lgm` -o SummersDiffCoeff

SummersDiffCoeffBounceAvg: SummersDiffCoeffBounceAvg.c DumpGif.c
	gcc SummersDiffCoeffBounceAvg.c DumpGif.c `pkg-config --cflags --libs lgm` -o SummersDiffCoeffBounceAvg

SummersDiffCoeffDerivBounceAvg: SummersDiffCoeffDerivBounceAvg.c DumpGif.c
	gcc SummersDiffCoeffDerivBounceAvg.c DumpGif.c `pkg-config --cflags --libs lgm` -o SummersDiffCoeffDerivBounceAvg


clean:
	rm SummersDiffCoeff SummersDiffCoeffBounceAvg SummersDiffCoeffDerivBounceAvg
