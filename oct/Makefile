OCTS = comp_wfac.oct comp_iwfac.oct comp_idgt_fac.oct			\
       comp_idgtreal_fac.oct comp_dgt_fb.oct comp_gabdual_long.oct	\
       comp_gabtight_long.oct comp_heapint.oct	\
       comp_col2diag.oct comp_gabreassign.oct comp_pgauss.oct		\
       comp_fftreal.oct comp_ifftreal.oct comp_dgt_long.oct		\
       comp_dwilt_long.oct comp_dgtreal_long.oct

all:	$(OCTS)

%.oct: %.cc config.h
	mkoctfile -Wall -I../thirdparty -I. -I../src -L../src -lltfat $<

clean:
	rm *.o *.oct

.PHONY: all clean