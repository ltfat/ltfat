CFLAGS=-O2 -Wall -std=c99 -I../src -L../lib
LIBS= -lltfat -lfftw3 -lm

timers = time_dgt_fac time_dgt_fb time_dgt_ola time_dgtreal_ola \
	time_dgtreal_fb time_dgtreal_fac time_fftw time_dgt_shear time_dgt_multi \
	time_dgt_baselinefft

all: 	$(timers)

time_%:	time_%.o ltfat_time.o
	$(CC) $(CFLAGS) $< ltfat_time.o $(LIBS) -o $@

%.o: %.c Makefile config.h
	$(CC) $(CFLAGS) -I../thirdparty -c $<

clean:
	rm *.o $(timers)
