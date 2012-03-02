# Very simple makefile illustrating how to use pkg-config to compile


all: ScatterMesh ScatterMesh2


ScatterMesh: ScatterMesh.c
	gcc ScatterMesh.c `pkg-config --cflags --libs lgm` -o ScatterMesh

ScatterMesh2: ScatterMesh2.c
	gcc ScatterMesh2.c `pkg-config --cflags --libs lgm` -o ScatterMesh2


clean:
	rm ScatterMesh ScatterMesh2
