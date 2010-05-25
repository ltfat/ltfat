CFLAGS=-O2 -Wall -std=c99 -I../src -L../lib
LIBS= -lltfat -lfftw3 -lm

all: time_dgt_fac time_dgt_fb

time_dgt_fb: time_dgt_fb.o ltfat_time.o
	$(CC) $(CFLAGS) time_dgt_fb.o ltfat_time.o $(LIBS) -o time_dgt_fb

time_dgt_fac: time_dgt_fac.o ltfat_time.o
	$(CC) $(CFLAGS) time_dgt_fac.o ltfat_time.o $(LIBS) -o time_dgt_fac

%.o: %.c Makefile config.h
	$(CC) $(CFLAGS) -I../thirdparty -c $<

clean:
	rm *.o time_dgt_fb time_dgt_fac
