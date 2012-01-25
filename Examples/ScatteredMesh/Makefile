# Very simple makefile illustrating how to use pkg-config to compile


ScatterMesh: ScatterMesh.c
	gcc ScatterMesh.c `pkg-config --cflags --libs lgm` -o ScatterMesh

clean:
	rm ScatterMesh
