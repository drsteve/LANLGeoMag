# Very simple makefile illustrating how to use pkg-config to compile

all: DynamicMemory

DynamicMemory: DynamicMemory.c
	gcc DynamicMemory.c `pkg-config --cflags --libs lgm` -Wall -o DynamicMemory

clean:
	rm DynamicMemory
