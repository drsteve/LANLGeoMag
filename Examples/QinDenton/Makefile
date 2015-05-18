# Very simple makefile illustrating how to use pkg-config to compile

all: QinDenton QinDenton2

QinDenton: QinDenton.c
	gcc QinDenton.c `pkg-config --cflags --libs lgm` -o QinDenton

QinDenton2: QinDenton2.c
	gcc QinDenton2.c `pkg-config --cflags --libs lgm` -o QinDenton2


clean:
	rm QinDenton QinDenton2
