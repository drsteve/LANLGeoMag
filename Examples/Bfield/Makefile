all: B GradB GradB2 B_Cross_GradB_Over_B CurlB CurlB2 DivB


B: B.c
	gcc B.c `pkg-config --cflags --libs lgm` -o B -Wall

GradB: GradB.c
	gcc GradB.c `pkg-config --cflags --libs lgm` -o GradB -Wall

GradB2: GradB2.c
	gcc GradB2.c `pkg-config --cflags --libs lgm` -o GradB2 -Wall

B_Cross_GradB_Over_B: B_Cross_GradB_Over_B.c
	gcc B_Cross_GradB_Over_B.c `pkg-config --cflags --libs lgm` -o B_Cross_GradB_Over_B -Wall

CurlB: CurlB.c
	gcc CurlB.c `pkg-config --cflags --libs lgm` -o CurlB -Wall

CurlB2: CurlB2.c
	gcc CurlB2.c `pkg-config --cflags --libs lgm` -o CurlB2 -Wall

DivB: DivB.c
	gcc DivB.c `pkg-config --cflags --libs lgm` -o DivB -Wall


clean:
	rm *~ *.o B GradB GradB2 B_Cross_GradB_Over_B CurlB CurlB2 DivB
