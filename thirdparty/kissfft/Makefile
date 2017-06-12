FILES = $(wildcard *.c)
OFILES = $(patsubst %.c,%.o,$(FILES))
DOFILES = $(addprefix d,$(OFILES))
SOFILES = $(addprefix s,$(OFILES))

all: $(DOFILES) $(SOFILES)

double: $(DOFILES)
single: $(SOFILES)
	
d%.o: %.c
	$(CC) -I../../include -c -DLTFAT_DOUBLE $< -o $@ 

s%.o: %.c
	$(CC) -I../../include -c -DLTFAT_SINGLE $< -o $@ 

.PHONY double single


