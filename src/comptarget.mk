ifeq ($(COMPTARGET),debug)
	# Do debug stuff here
else
	CFLAGS += -O2 -DNDEBUG
endif
