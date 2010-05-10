# Use GNU Make to process this file
CC = gcc
LCCLIB = '/c/Program Files/MATLAB/R2009b/sys/lcc/bin/lcclib'

CFLAGS=-O2 -Wall -fPIC -std=c99 -I../thirdparty -L../thirdparty

files = sdgt.o sdgt_fac.o sdgtreal_fac.o sidgt_fac.o sdgt_fb.o		\
	sdgt_walnut.o ddgt.o ddgt_fac.o ddgtreal_fac.o didgt_fac.o	\
	ddgt_fb.o ddgt_walnut.o sgabdual_fac.o sgabtight_fac.o		\
	sltfat_blas.o dgabdual_fac.o dgabtight_fac.o		\
	dltfat_blas.o sltfat_lapack.o swfac.o siwfac.o	\
	sreassign.o swindows.o dltfat_lapack.o dwfac.o diwfac.o		\
	dreassign.o dwindows.o sspread.o sheapint.o swinmanip.o		\
	dspread.o dheapint.o dwinmanip.o \
	sdwilt.o sgabdual.o sgabtight.o stfutil.o \
	ddwilt.o dgabdual.o dgabtight.o dtfutil.o \
	integer_manip.o 

all: libltfat.a unixnomem

winnomem: $(files)
	$(LCCLIB) /out:libltfat-nomem.lib $(files)
	mv -f libltfat-nomem.lib ../lib

unixnomem: $(files) 
	ar rvu libltfat-nomem.a $(files)
	ranlib libltfat-nomem.a
	cp -f libltfat-nomem.a ../lib

libltfat.a: $(files) c-safe-memalloc.o
	ar rvu libltfat.a $(files) c-safe-memalloc.o
	ranlib libltfat.a
	cp -f libltfat.a ../lib

s%.o: %.c config.h
	$(CC) $(CFLAGS) -DLTFAT_SINGLE  -c $< -o s$*.o

d%.o: %.c config.h
	$(CC) $(CFLAGS) -DLTFAT_DOUBLE  -c $< -o d$*.o

%.o: %.c Makefile config.h
	$(CC) $(CFLAGS) -I../thirdparty -c $<

clean:
	rm *.o *.a 
