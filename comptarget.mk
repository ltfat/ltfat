ifeq ($(COMPTARGET),debug)
	# Do debug stuff here
	CFLAGS += -O0 -g
	CXXFLAGS += -O0 -g
else
	CFLAGS += -O2 -DNDEBUG
	CXXFLAGS += -O2 -DNDEBUG
endif
