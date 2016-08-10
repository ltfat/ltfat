
CFILES = $(shell ls test_*.c)
HFILES = $(shell ls test_*.c)

test_all_libltfat: Makefile ../../build/libltfat.a $(HFILES) $(CFILES)
	$(CC) -Wall -Wextra -pedantic -std=gnu99 -O0 -g -I../../include -I../../thirdparty test_all_libltfat.c -o test_all_libltfat ../../build/libltfat.a -lfftw3 -lfftw3f -lm

run: test_all_libltfat
	./test_all_libltfat

mem: test_all_libltfat
	valgrind --leak-check=yes ./test_all_libltfat 

clean:
	-rm -f test_all_libltfat

