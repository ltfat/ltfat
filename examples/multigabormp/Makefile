CFLAGS+=-Ofast -Wall -Wextra -std=c++11 -DLTFAT_SINGLE

SRC=$(wildcard *.cpp)
PROGS = $(patsubst %.cpp,%,$(SRC))

all: $(PROGS)

%: %.cpp
	$(CXX) $(CFLAGS) -I.. $< -o $@ ../../build/libltfat.a -lblas -llapack -lfftw3 -lfftw3f -lc -lm -lsndfile

clean:
	-rm $(PROGS)

