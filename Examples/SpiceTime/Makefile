# Very simple makefile illustrating how to use pkg-config to compile


SpiceTime: SpiceTime.c
	gcc SpiceTime.c `pkg-config --cflags --libs lgm` -o SpiceTime

clean:
	rm SpiceTime
