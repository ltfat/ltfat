# This is the main Makefile for libltfat
#
# Builds three static and three shared libraries (default prefix is build/):
#
# 	libltfat.a     Contains double, single and common code
# 	libltfatd.a    Contains double and common code
# 	libltfatf.a    Contains single and common code
#
# make CROSS=x86_64-w64-mingw32.static-
# or
# make CROSS=x86_64-w64-mingw32.static- NOBLASLAPACK=1
#

CC=$(CROSS)gcc
AR=$(CROSS)ar
OBJCOPY=$(CROSS)objcopy

# Base CFLAGS
CFLAGS+=-Ithirdparty -Wall -Wextra -pedantic -std=gnu99 $(OPTCFLAGS)

# Dependencies
LFLAGSBASE?=
ifndef NOBLASLAPACK
	LFLAGSBASE +=-llapack -lblas
endif

# The following adds parameters to CFLAGS
include comptarget.mk

# Define source files from src/
include filedefs.mk

ifdef CROSS
	buildprefix ?= build/$(CROSS)
	objprefix ?= obj/$(CROSS)
	EXTRALFLAGS = -Wl,--out-implib,$@.a -static-libgcc
	LFLAGSBASE +=-lgfortran -lquadmath
else
	CFLAGS += -fPIC
	buildprefix ?= build
	objprefix ?= obj
	EXTRALFLAGS =
endif

DLFLAGS = -Wl,--no-undefined $(LFLAGSBASE) -lfftw3 -lm
SLFLAGS =  -Wl,--no-undefined $(LFLAGSBASE) -lfftw3f -lm
ALLLFLAGS = -Wl,--no-undefined $(LFLAGSBASE) -lfftw3 -lfftw3f -lm
# Convert *.c names to *.o
toCompile = $(patsubst %.c,%.o,$(files))
toCompile_complextransp = $(patsubst %.c,%.o,$(files_complextransp))
toCompile_notypechange = $(patsubst  %.c,%.o,$(files_notypechange))

ifndef NOBLASLAPACK
	toCompile += $(patsubst %.c,%.o,$(files_blaslapack))
	toCompile_complextransp += $(patsubst %.c,%.o,$(files_blaslapack_complextransp))
endif

COMMONFILES = $(addprefix $(objprefix)/common/d,$(toCompile_notypechange))
COMMONFILESFORSFILES = $(addprefix $(objprefix)/common/s,$(toCompile_notypechange))

DFILES   = $(addprefix $(objprefix)/double/,$(toCompile) $(toCompile_complextransp)) \
		   $(addprefix $(objprefix)/complexdouble/,$(toCompile_complextransp))
SFILES   = $(addprefix $(objprefix)/single/,$(toCompile) $(toCompile_complextransp)) \
		   $(addprefix $(objprefix)/complexsingle/,$(toCompile_complextransp))

# Define targets
DTARGET=$(buildprefix)/libltfatd.a
STARGET=$(buildprefix)/libltfatf.a
ALLTARGET=$(buildprefix)/libltfat.a

ifndef CROSS
	SO_DTARGET=$(patsubst %.a,%.so,$(DTARGET))
	SO_STARGET=$(patsubst %.a,%.so,$(STARGET))
	SO_ALLTARGET=$(patsubst %.a,%.so,$(ALLTARGET))
else
	SO_DTARGET=$(patsubst %.a,%.dll,$(DTARGET))
	SO_STARGET=$(patsubst %.a,%.dll,$(STARGET))
	SO_ALLTARGET=$(patsubst %.a,%.dll,$(ALLTARGET))
endif

all: $(DTARGET) $(STARGET) $(SO_DTARGET) $(SO_STARGET) $(ALLTARGET) $(SO_ALLTARGET)

$(ALLTARGET): $(STARGET) $(DTARGET)
	ar rvu $@ $(COMMONFILES) $(DFILES) $(SFILES)
	ranlib $@

$(DTARGET): $(buildprefix) $(objprefix)/double $(objprefix)/complexdouble $(objprefix)/common $(DFILES) $(COMMONFILES)
	ar rvu $@ $(DFILES) $(COMMONFILES)
	ranlib $@

$(STARGET): $(buildprefix) $(objprefix)/single $(objprefix)/complexsingle $(objprefix)/common $(SFILES) $(COMMONFILESFORSFILES)
	ar rvu $@ $(SFILES) $(COMMONFILESFORSFILES)
	ranlib $@

$(SO_ALLTARGET): $(ALLTARGET)
	$(CC) -shared -fPIC -o $@ $(COMMONFILES) $(DFILES) $(SFILES)  $(ALLLFLAGS) $(EXTRALFLAGS)

$(SO_DTARGET): $(DTARGET)
	$(CC) -shared -fPIC -o $@ $(COMMONFILES) $(DFILES) $(DLFLAGS) $(EXTRALFLAGS)

$(SO_STARGET): $(STARGET)
	$(CC) -shared -fPIC -o $@ $(COMMONFILESFORSFILES) $(SFILES) $(SLFLAGS) $(EXTRALFLAGS)

$(objprefix)/common/d%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@ 

$(objprefix)/common/s%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@ 
	# Overwrite symbols to avoid dependency on fftw since we actually link with fftwf
	$(OBJCOPY) --redefine-sym fftw_malloc=fftwf_malloc $@
	$(OBJCOPY) --redefine-sym fftw_free=fftwf_free $@

$(objprefix)/common/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(OPTCFLAGS)

$(objprefix)/single/%.o: src/%.c
	$(CC) $(CFLAGS) -DLTFAT_SINGLE  -c $< -o $@

$(objprefix)/double/%.o: src/%.c
	$(CC) $(CFLAGS) -DLTFAT_DOUBLE  -c $< -o $@

$(objprefix)/complexsingle/%.o: src/%.c
	$(CC) $(CFLAGS) -DLTFAT_SINGLE -DLTFAT_COMPLEXTYPE -c $< -o $@

$(objprefix)/complexdouble/%.o: src/%.c
	$(CC) $(CFLAGS) -DLTFAT_DOUBLE -DLTFAT_COMPLEXTYPE -c $< -o $@

$(buildprefix):
	@mkdir -p $(buildprefix)

$(objprefix)/common:
	@mkdir -p $(objprefix)/common

$(objprefix)/double:
	@mkdir -p $(objprefix)/double

$(objprefix)/single:
	@mkdir -p $(objprefix)/single

$(objprefix)/complexdouble:
	@mkdir -p $(objprefix)/complexdouble

$(objprefix)/complexsingle:
	@mkdir -p $(objprefix)/complexsingle

.PHONY: clean help

clean:
	@rm -rf build obj

help:
	@echo "USAGE: make [target]"
	@echo "Options:"
	@echo "    make [target] CONFIG=debug               Compiles the library in a debug mode"
	@echo "    make [target] NOBLASLAPACK=              Compiles the library without BLAS and LAPACK dependencies"
