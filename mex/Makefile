# Use GNU Make to process this file.

EXT = $(shell mexext)
# This line gets the MEX extension for the platform. If the shell
# command 'mexext' is not available on your Matlab installation,
# please execute the 'mexext' command within Matlab, and put the
# answer into the line below (remember to uncomment).
#EXT = $(glnx86)

MEXS = comp_wfac comp_iwfac comp_dgt_long comp_idgt_fac			\
       comp_idgtreal_fac comp_dgt_fb comp_gabdual_long			\
       comp_gabtight_long comp_heapint comp_col2diag comp_gabreassign	\
       comp_pgauss comp_fftreal comp_ifftreal comp_dwilt_long

MEXWITHEXTS = $(addsuffix .$(EXT),$(MEXS))

CFLAGS=-I../thirdparty -I. -I../src -L../lib


all:	$(MEXWITHEXTS)

%.$(EXT): %.c Makefile config.h ../src/libltfat.a ../src/ltfat.h
	mex $(CFLAGS) $< -lfftw3 -lfftw3f -lblas -llapack -lltfat

clean:
	rm *.$(EXT)

.PHONY: all clean
