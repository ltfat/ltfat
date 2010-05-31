# Use GNU Make to process this file
CC = gcc
LCCLIB = '/c/Program Files/MATLAB/R2009b/sys/lcc/bin/lcclib'

CFLAGS=-O2 -Wall -fPIC -std=c99 -I../thirdparty -L../thirdparty

files = sdgt.o sdgt_fac.o sdgtreal_fac.o sidgt_fac.o sdgt_fb.o		\
	sdgt_walnut.o ddgt.o ddgt_fac.o ddgtreal_fac.o didgt_fac.o	\
	ddgt_fb.o ddgt_walnut.o sgabdual_fac.o sgabtight_fac.o		\
	dgabdual_fac.o dgabtight_fac.o		\
	swfac.o siwfac.o	\
	sreassign.o swindows.o dwfac.o diwfac.o		\
	dreassign.o dwindows.o sspread.o sheapint.o swinmanip.o		\
	dspread.o dheapint.o dwinmanip.o \
	sdwilt.o sgabdual.o sgabtight.o stfutil.o \
	ddwilt.o dgabdual.o dgabtight.o dtfutil.o \
	integer_manip.o 

files_unix = $(files) dltfat_blaslapack.o sltfat_blaslapack.o
files_matlab = $(files) dltfat_blaslapack_matlab.o sltfat_blaslapack_matlab.o

all: libltfat.a unixnomem

winnomem: $(files_matlab) 
	$(LCCLIB) /out:libltfat-nomem.lib $(files_matlab)
	mv -f libltfat-nomem.lib ../lib

unixnomem: $(files_unix) 
	ar rvu libltfat-nomem.a $(files_unix)
	ranlib libltfat-nomem.a
	cp -f libltfat-nomem.a ../lib

libltfat.a: $(files_unix) c-safe-memalloc.o 
	ar rvu libltfat.a $(files_unix) c-safe-memalloc.o
	ranlib libltfat.a
	cp -f libltfat.a ../lib

sltfat_blaslapack_matlab.o: ltfat_blaslapack.c config.h
	$(CC) $(CFLAGS) -DLTFAT_SINGLE -DMATLABFORTRAN -c $< -o $*.o

dltfat_blaslapack_matlab.o: ltfat_blaslapack.c config.h
	$(CC) $(CFLAGS) -DLTFAT_DOUBLE -DMATLABFORTRAN -c $< -o $*.o


s%.o: %.c config.h
	$(CC) $(CFLAGS) -DLTFAT_SINGLE  -c $< -o s$*.o

d%.o: %.c config.h
	$(CC) $(CFLAGS) -DLTFAT_DOUBLE  -c $< -o d$*.o

%.o: %.c Makefile config.h
	$(CC) $(CFLAGS) -I../thirdparty -c $<

clean:
	rm *.o *.a 
