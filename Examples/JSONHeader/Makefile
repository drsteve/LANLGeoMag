all: writeJSONHeader

writeJSONHeader: writeJSONHeader.c
	gcc writeJSONHeader.c `pkg-config --cflags --libs lgm` -g -o writeJSONHeader -Wall


clean:
	rm -f *~ *.o writeJSONHeader
