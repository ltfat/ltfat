

FILES = $(wildcard *.c)
OFILES = $(patsubst %.c,%.o,$(FILES))


all: $(OFILES)
	
%.o: %.c
	$(CC) -c -Dkiss_fft_scalar=double $< -o $@ 


