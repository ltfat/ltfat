CFLAGS+=-g -Wall -Wextra -std=gnu99 -O3

SRC=$(wildcard *.c)
PROGS = $(patsubst %.c,%,$(SRC))

all: $(PROGS)

%: %.c
	$(CC) $(CFLAGS) -g -Ofast -I.. $< -o $@ ../build/libltfat.a -lblas -llapack -lfftw3 -lc -lm -lsndfile

clean:
	-rm $(PROGS)
