# Very simple makefile illustrating how to use pkg-config to compile

CC=gcc
PKG_CONFIG=pkg-config

quartic: quartic.c Makefile
	$(CC) quartic.c `$(PKG_CONFIG) --cflags --libs lgm` -Wall -o quartic

clean:
	rm quartic
