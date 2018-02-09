CFLAGS+=-g -Wall -Wextra -std=c++11
CPP=g++

SRC=$(wildcard *.c)
PROGS = $(patsubst %.c,%,$(SRC))

all: $(PROGS)

%: %.c
	$(CPP) $(CFLAGS) -g -Ofast -I.. $< -o $@ ../build/libltfat.a -lblas -llapack -lfftw3 -lc -lm -lsndfile

clean:
	-rm $(PROGS)
