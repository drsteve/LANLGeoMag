all: writeJSONHeader

writeJSONHeader: writeJSONHeader.c
	gcc writeJSONHeader.c `pkg-config --cflags --libs lgm` -o writeJSONHeader -Wall


clean:
	rm *~ *.o writeJSONHeader
