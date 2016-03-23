# Very simple makefile illustrating how to use pkg-config to compile
LIBS   = `pkg-config --cflags --libs lgm`

all: JulianDate LeapSecondDateBoundary TestLeapSeconds TimeConv TimeConv2 TimeConv3 TimeConvOnLeapSecondDay

JulianDate: JulianDate.c
	gcc JulianDate.c $(LIBS) -o JulianDate

LeapSecondDateBoundary: LeapSecondDateBoundary.c
	gcc LeapSecondDateBoundary.c $(LIBS) -o LeapSecondDateBoundary

TestLeapSeconds: TestLeapSeconds.c
	gcc TestLeapSeconds.c $(LIBS) -o TestLeapSeconds

TimeConv: TimeConv.c
	gcc TimeConv.c $(LIBS) -o TimeConv

TimeConv2: TimeConv2.c
	gcc TimeConv2.c $(LIBS) -o TimeConv2

TimeConv3: TimeConv3.c
	gcc TimeConv3.c $(LIBS) -o TimeConv3

TimeConvOnLeapSecondDay: TimeConvOnLeapSecondDay.c
	gcc TimeConvOnLeapSecondDay.c $(LIBS) -o TimeConvOnLeapSecondDay

clean:
	rm JulianDate LeapSecondDateBoundary TestLeapSeconds TimeConv TimeConv2 TimeConv3 TimeConvOnLeapSecondDay 
