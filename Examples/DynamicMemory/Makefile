# Very simple makefile illustrating how to use pkg-config to compile

CC=gcc
PKG_CONFIG=pkg-config

all: DynamicMemory

DynamicMemory: DynamicMemory.c
	$(CC) DynamicMemory.c `$(PKG_CONFIG) --cflags --libs lgm` -Wall -o DynamicMemory

clean:
	rm DynamicMemory
