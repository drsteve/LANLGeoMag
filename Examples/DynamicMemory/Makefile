# Very simple makefile illustrating how to use pkg-config to compile

CC=gcc
PKG_CONFIG=pkg-config

all: DynamicMemory

DynamicMemory: DynamicMemory.c
	$(CC) DynamicMemory.c `$(PKG_CONFIG) --cflags --libs lgm` `$(PKG_CONFIG) --cflags --libs hdf5 &2> /dev/null` -Wall -o DynamicMemory

clean:
	rm DynamicMemory
