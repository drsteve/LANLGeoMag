LGMFLAGS = `pkg-config lgm --cflags --libs`


all: ZeroTilt


ZeroTilt: ZeroTilt.c
	h5cc ZeroTilt.c $(LGMFLAGS) -o ZeroTilt -Wall

clean:
	rm *~ *.o ZeroTilt
