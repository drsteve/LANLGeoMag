all: UBK


Vec: UBK.c
	gcc UBK.c `pkg-config --cflags --libs lgm` -o UBK -Wall



clean:
	rm -f *~ *.o UBK
