CFLAGS+=-g -Wall -Wextra -std=c++11 
CPP=g++

SRC=$(wildcard *.cpp)
PROGS = $(patsubst %.cpp,%,$(SRC))

all: $(PROGS)

%: %.c
	$(CPP) $(CFLAGS) -g -Ofast -I.. -I../../include $< -o $@ ../../build/libltfat.a -lblas -llapack -lfftw3 -lc -lm -lsndfile

clean:
	-rm $(PROGS)
