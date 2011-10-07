all: Vcg Vcg2

Vcg: Vcg.c
	gcc Vcg.c `pkg-config --cflags --libs lgm` -o Vcg -Wall

Vcg2: Vcg2.c
	gcc Vcg2.c `pkg-config --cflags --libs lgm` -o Vcg2 -Wall


clean:
	rm *~ *.o Vcg Vcg2
