HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`


all: B Bcompare B2 GradBvec GradBvec2 GradB GradB2 B_Cross_GradB_Over_B CurlB CurlB2 DivB


B: B.c
	gcc B.c $(LGMFLAGS) $(HDF5FLAGS) -o B -Wall

Bcompare: Bcompare.c
	gcc Bcompare.c $(LGMFLAGS) $(HDF5FLAGS) -o Bcompare -Wall

B2: B2.c
	gcc B2.c $(LGMFLAGS) $(HDF5FLAGS) -o B2 -Wall

GradBvec: GradBvec.c
	gcc GradBvec.c $(LGMFLAGS) $(HDF5FLAGS) -o GradBvec -Wall

GradBvec2: GradBvec2.c
	gcc GradBvec2.c $(LGMFLAGS) $(HDF5FLAGS) -o GradBvec2 -Wall

GradB: GradB.c
	gcc GradB.c $(LGMFLAGS) $(HDF5FLAGS) -o GradB -Wall

GradB2: GradB2.c
	gcc GradB2.c $(LGMFLAGS) $(HDF5FLAGS) -o GradB2 -Wall

B_Cross_GradB_Over_B: B_Cross_GradB_Over_B.c
	gcc B_Cross_GradB_Over_B.c $(LGMFLAGS) $(HDF5FLAGS) -o B_Cross_GradB_Over_B -Wall

CurlB: CurlB.c
	gcc CurlB.c $(LGMFLAGS) $(HDF5FLAGS) -o CurlB -Wall

CurlB2: CurlB2.c
	gcc CurlB2.c $(LGMFLAGS) $(HDF5FLAGS) -o CurlB2 -Wall

DivB: DivB.c
	gcc DivB.c $(LGMFLAGS) $(HDF5FLAGS) -o DivB -Wall


clean:
	rm *~ *.o B Bcompare B2 GradBvec GradBvec2 GradB GradB2 B_Cross_GradB_Over_B CurlB CurlB2 DivB
