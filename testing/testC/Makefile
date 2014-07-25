CFLAGS=-g -Wall -Wextra -DNDEBUG
LIBS=

TEST_SRC=$(wildcard test_*.c)
TESTS=$(patsubst %.c,%,$(TEST_SRC))

all: $(TESTS)

$(TESTS): % : %.c
	$(CC) $(CFLAGS) $^ $(shell sed '2q;d' $^) -o $@


debug: CFLAGS=-g -Wall -Wextra
debug: all



clean:
	rm -rf $(TESTS)

.PHONY: tests
tests: $(TESTS)
	sh ./runtests.sh

valgrind:
	VALGRIND="valgrind --log-file=/tmp/valgrind-%p.log" $(MAKE)
