all: B B2


B: B.c
	gcc B.c `pkg-config --cflags --libs lgm` -o B -Wall

B2: B2.c
	gcc B2.c `pkg-config --cflags --libs lgm` -o B2 -Wall

clean:
	rm B B2
