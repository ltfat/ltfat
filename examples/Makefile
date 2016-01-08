CFLAGS+=-g -Wall -Wextra -std=gnu99

testgab: example_gabanasyn.c
	$(CC) $(CFLAGS) -I.. example_gabanasyn.c -o example_gabanasyn ../libltfat.a -lblas -llapack -lfftw3 -lc -lm 
