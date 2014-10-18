# Very simple makefile illustrating how to use pkg-config to compile

all: LstarVersusPA LCDS

LstarVersusPA : LstarVersusPA.c
	gcc LstarVersusPA.c `pkg-config --cflags --libs lgm` -o LstarVersusPA

LCDS: LCDS.c
	gcc LCDS.c `pkg-config --cflags --libs lgm` -o LCDS

clean:
	rm LCDS LstarVersusPA
