# Very simple makefile illustrating how to use pkg-config to compile

all: FluxTubeVolume

FluxTubeVolume: FluxTubeVolume.c
	gcc FluxTubeVolume.c `pkg-config --cflags --libs lgm` -o FluxTubeVolume

clean:
	rm FluxTubeVolume
