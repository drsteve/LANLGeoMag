all: Vec


Vec: Vec.c
	gcc Vec.c `pkg-config --cflags --libs lgm` -o Vec -Wall



clean:
	rm -f *~ *.o Vec