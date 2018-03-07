CFLAGS+=-Ofast -Wall -Wextra -std=c++1z -DLTFAT_DOUBLE -DCXXOPTS_USE_UNICODE

SRC=$(wildcard *.cpp)
PROGS = $(patsubst %.cpp,%,$(SRC))

all: $(PROGS)

%: %.cpp
	$(CXX) $(CFLAGS) -I../utils $< -o $@ ../../build/libltfat.a -lblas -llapack -lfftw3 -lfftw3f -lc -lm -lsndfile -licuuc

clean:
	-rm $(PROGS)

