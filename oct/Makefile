# Use GNU Make to process this file.

OCTSCC   = $(shell ls comp_*.cc)
OCTSBASE = $(basename $(OCTSCC))
OCTS     = $(addsuffix .oct,$(OCTSBASE))

all:	$(OCTS)

%.oct: %.cc config.h
	mkoctfile -Wall -I../thirdparty -I. -I../src -L../src -lltfat $<

clean:
	rm *.o *.oct

.PHONY: all clean