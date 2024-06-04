#CFLAGS    = -Wall -O3 -g -c -Wall -mfpmath=387 -funroll-all-loops -ffast-math -fopenmp
#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
#FFLAGS    =  -fno-automatic 
LIBS      = `pkg-config --cflags --libs lgm`
CC        = gcc

FILES1    =  

all   : TwoLineElements TwoLineElementsFromFile NearestTLE Forte_Meteor

Forte_Meteor : Makefile $(FILES1) Forte_Meteor.c
	$(CC) $(FILES1) Forte_Meteor.c $(LIBS) -o Forte_Meteor

TwoLineElements : Makefile $(FILES1) TwoLineElements.c
	$(CC) $(FILES1) TwoLineElements.c $(LIBS) -o TwoLineElements

TwoLineElementsFromFile : Makefile $(FILES1) TwoLineElementsFromFile.c
	$(CC) $(FILES1) TwoLineElementsFromFile.c $(LIBS) -o TwoLineElementsFromFile

NearestTLE : Makefile $(FILES1) NearestTLE.c
	$(CC) $(FILES1) NearestTLE.c $(LIBS) -o NearestTLE


clean :
	rm TwoLineElements TwoLineElementsFromFile NearestTLE *.o
