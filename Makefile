# Use GNU Make to process this file
#CC = i686-pc-mingw32-gcc 
CC = gcc


CFLAGS=-O2 -Wall -fPIC -std=c99 -I../thirdparty -L../thirdparty

files = \
	sdgt.o sdgt_fac.o sdgtreal_fac.o sidgt_fac.o sdgt_fb.o		\
	sdgt_walnut.o ddgt.o ddgt_fac.o ddgtreal_fac.o didgt_fac.o	\
	ddgt_ola.o sdgt_ola.o \
	ddgt_fb.o ddgt_walnut.o \
	swfac.o siwfac.o sidgt_fb.o didgt_fb.o \
	sreassign.o swindows.o dwfac.o diwfac.o		\
	dreassign.o dwindows.o sspread.o sheapint.o swinmanip.o		\
	dspread.o dheapint.o dwinmanip.o \
	sdwilt.o stfutil.o sfilterbank.o \
	ddwilt.o dtfutil.o dfilterbank.o \
	spfilt.o dpfilt.o \
	integer_manip.o 

files_blaslapack = \
	sgabdual.o sgabtight.o sgabdual_fac.o sgabtight_fac.o \
	dgabdual.o dgabtight.o dgabdual_fac.o dgabtight_fac.o 	

files_unix = $(files) $(files_blaslapack) dltfat_blaslapack.o sltfat_blaslapack.o
files_matlab = $(files) $(files_blaslapack) dltfat_blaslapack_matlab.o sltfat_blaslapack_matlab.o

all: libltfat.a unixnomem

win_ms: $(files) c-safe-memalloc.o
	gcc -shared -o ltfat.dll -Wl,--output-def,ltfat.def,--out-implib,libltfat_dll.a \
		-L../thirdparty -lfftw3-3 -lfftw3f-3 $(files) c-safe-memalloc.o
	cp -f ltfat.dll ../mex/
	cp -f ltfat.def ../mex/
	echo "To finish the creation, please run the Microsoft lib tool on ltfat.dll and ltfat.def"

wincrosscompile: $(files_matlab)
	ar rvu libltfat-nomem.lib $(files_matlab)
	ranlib libltfat-nomem.lib
	cp -f libltfat-nomem.lib ../lib/

winnomem: $(files_matlab)
	ar rvu libltfat-nomem.lib $(files_matlab)
	ranlib libltfat-nomem.lib
	dlltool -z libltfat-nomem.def --export-all-symbols libltfat-nomem.lib
	cp -f libltfat-nomem.lib ../lib/
	cp -f libltfat-nomem.def ../lib/

unixnomem: $(files_unix) 
	ar rvu libltfat-nomem.a $(files_unix)
	ranlib libltfat-nomem.a
	cp -f libltfat-nomem.a ../lib

libltfat.a: $(files_unix) c-safe-memalloc.o 
	ar rvu libltfat.a $(files_unix) c-safe-memalloc.o
	ranlib libltfat.a
	cp -f libltfat.a ../lib

libltfat.so: $(files_unix) c-safe-memalloc.o
	gcc -shared -Wl,-soname,libltfat.so.1 -o libltfat.so.1.0 $(files_unix) c-safe-memalloc.o

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
