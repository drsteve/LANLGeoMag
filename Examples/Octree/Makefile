# Very simple makefile illustrating how to use pkg-config to compile


kNN: kNN.c
	gcc kNN.c `pkg-config --cflags --libs lgm` -o kNN

clean:
	rm kNN
