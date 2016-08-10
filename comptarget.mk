ifeq ($(COMPTARGET),debug)
	# Do debug stuff here
	CFLAGS += -O0 -g
	CPPFLAGS += -O0 -g
else
	CFLAGS += -O2 -DNDEBUG
	CPPFLAGS += -O2 -DNDEBUG
endif
