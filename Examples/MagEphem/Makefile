#CFLAGS    = -Wall -O3 -march=opteron -msse2 -mfpmath=sse,387 -ffast-math
CDF_LIBS   = -lcdf -L/usr/local/cdf/lib
CDF_INC    = -I/usr/local/cdf/include

LIBS      := `pkg-config --libs lgm` $(CDF_LIBS) `pkg-config --libs hdf5 2> /dev/null`
CFLAGS    := `pkg-config --cflags lgm` $(CDF_INC) -g `pkg-config --cflags hdf5 2> /dev/null`
CC        = gcc -Wno-unused-variable
FILES     = MagEphemFromFile

all: $(FILES)

-include extra_flags

extra_flags:
	if [[ "${MACHTYPE}" =~ "*darwin*" ]]; then echo "LIBS := -largp -I/opt/local/include -L/opt/local/lib" > extra_flags; echo "CFLAGS := -I/opt/local/include -L/opt/local/lib -largp" >> extra_flags; else echo > extra_flags; fi

$(FILES): %: %.c Makefile
	$(CC) $< $(LIBS) $(CFLAGS) -o $@

clean :
	rm MagEphemFromFile *.o
