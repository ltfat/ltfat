CFLAGS=-Wall

all:	mex_dgt_fac_1.oct mex_dgt_fac_2.oct mex_dgt_fac_3.oct \
	mex_dgt_fac_4.oct mex_dgt_fac_5.oct mex_dgt_fac_6.oct \
	mex_dgt_fac_7.oct mex_dgt_fb_1.oct mex_dgt_fb_2.oct

%.oct: %.cc Makefile
	mkoctfile $(CFLAGS) -I. -I.. $<

clean:	
	rm *.o *.oct

.PHONY: all clean
