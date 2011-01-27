# Very simple makefile illustrating how to use pkg-config to compile

all: CgmCoords

CgmCoords: CgmCoords.c
	gcc CgmCoords.c `pkg-config --cflags --libs lgm` -o CgmCoords


clean:
	rm CgmCoords
