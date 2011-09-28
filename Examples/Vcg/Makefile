all: Vcg

Vcg: Vcg.c
	gcc Vcg.c `pkg-config --cflags --libs lgm` -o Vcg -Wall

clean:
	rm *~ *.o Vcg
