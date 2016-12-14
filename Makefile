# This is the main Makefile for libltfat
#
# Builds three static and three shared libraries (default prefix is build/):
#
# 	libltfat.a(.so)     Contains double, single and common code
# 	libltfatd.a(.so)    Contains double and common code
# 	libltfatf.a(.so)    Contains single and common code
#
# make CROSS=x86_64-w64-mingw32.static-
# or
# make CROSS=x86_64-w64-mingw32.static- NOBLASLAPACK=1
#
# Examples:
# ---------
#
# Linking against BLAS and LAPACK from Matlab (if shared libs are desired):
#  		make BLASLAPACKLIBS="-L/usr/local/MATLAB_R2013a/bin/glnxa64 -lmwblas -lmwlapack" 
#
#  		make clean && make BLASLAPACKLIBS="-L/usr/local/MATLAB_R2013a/bin/glnxa64 -lmwblas -lmwlapack" && make build/ltfat.h

include ostools.mk

ifdef CROSS
	CC=$(CROSS)gcc
	CXX=$(CROSS)g++
	AR=$(CROSS)ar
	OBJCOPY=$(CROSS)objcopy
	RANLIB=$(CROSS)ranlib
	buildprefix ?= build/$(CROSS)
	objprefix ?= obj/$(CROSS)
	MINGW=1
else
	CC?=gcc
	CXX?=g++
	AR?=ar
	OBJCOPY?=objcopy
	RANLIB?=ranlib
	buildprefix ?= build
	objprefix ?= obj
endif

MATLABPATH ?= /usr/local/MATLAB_R2013a/bin/glnxa64

PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

# Base CFLAGS
CFLAGS+=-Wall -Wextra -pedantic -std=c99 -Iinclude -Ithirdparty $(OPTCFLAGS)
CXXFLAGS+=-Wall -Wextra -pedantic -std=c++11 -fno-exceptions -fno-rtti -Iinclude -Ithirdparty $(OPTCFLAGS)

# The following adds parameters to CFLAGS
COMPTARGET ?= release
include comptarget.mk

# Define source files from src/
include filedefs.mk

FFTWLIBS?=-lfftw3 -lfftw3f
LFLAGS = -Wl,--no-undefined $(OPTLPATH)

ifdef MINGW
	EXTRALFLAGS = -Wl,--out-implib,$@.a -static-libgcc
	BLASLAPACKLIBS?=-llapack -lblas -lgfortran -lquadmath
	CFLAGS += -DLTFAT_BUILD_SHARED
	CXXFLAGS += -DLTFAT_BUILD_SHARED
else
	CFLAGS += -fPIC
	CXXFLAGS += -fPIC
	BLASLAPACKLIBS?=-llapack -lblas
endif

ifdef USECPP
ifeq ($(USECPP),1)
	CC = $(CXX)
	CFLAGS = $(CXXFLAGS)
endif
endif

# Dependencies
ifndef NOBLASLAPACK
	LFLAGS += $(BLASLAPACKLIBS)
endif
LFLAGS += $(FFTWLIBS) $(EXTRALFLAGS) $(OPTLFLAGS) -lm

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

# Define libraries
DSTATIC = libltfatd.a
SSTATIC = libltfatf.a
DSSTATIC = libltfat.a

ifndef MINGW
	DSHARED = $(patsubst %.a,%.so,$(DSTATIC))
	SSHARED = $(patsubst %.a,%.so,$(SSTATIC))
	DSSHARED = $(patsubst %.a,%.so,$(DSSTATIC))
else
	DSHARED = $(patsubst %.a,%.dll,$(DSTATIC))
	SSHARED = $(patsubst %.a,%.dll,$(SSTATIC))
	DSSHARED = $(patsubst %.a,%.dll,$(DSSTATIC))
endif

# Define targets
DTARGET=$(buildprefix)/$(DSTATIC)
STARGET=$(buildprefix)/$(SSTATIC)
DSTARGET=$(buildprefix)/$(DSSTATIC)
SO_DTARGET=$(buildprefix)/$(DSHARED)
SO_STARGET=$(buildprefix)/$(SSHARED)
SO_DSTARGET=$(buildprefix)/$(DSSHARED)

DDEP = $(buildprefix) $(objprefix)/double $(objprefix)/complexdouble $(objprefix)/common
SDEP = $(buildprefix) $(objprefix)/single $(objprefix)/complexsingle $(objprefix)/common

all: static shared

$(DSTARGET): $(DDEP) $(SDEP) $(COMMONFILES) $(DFILES) $(SFILES)
	$(AR) rv $@ $(COMMONFILES) $(DFILES) $(SFILES)
	$(RANLIB) $@

$(DTARGET): $(DDEP) $(DFILES) $(COMMONFILES)
	$(AR) rv $@ $(DFILES) $(COMMONFILES)
	$(RANLIB) $@

$(STARGET): $(SDEP) $(SFILES) $(COMMONFILESFORSFILES)
	$(AR) rv $@ $(SFILES) $(COMMONFILESFORSFILES)
	$(RANLIB) $@

$(SO_DSTARGET): $(DDEP) $(SDEP) $(COMMONFILES) $(DFILES) $(SFILES)
	$(CC) -shared -fPIC -o $@ $(COMMONFILES) $(DFILES) $(SFILES) $(LFLAGS)

$(SO_DTARGET): $(DDEP) $(COMMONFILES) $(DFILES)
	$(CC) -shared -fPIC -o $@ $(COMMONFILES) $(DFILES) $(LFLAGS)

$(SO_STARGET): $(SDEP) $(SFILES) $(COMMONFILESFORSFILES)
	$(CC) -shared -fPIC -o $@ $(COMMONFILESFORSFILES) $(SFILES) $(LFLAGS)

$(objprefix)/common/d%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_DOUBLE -c $< -o $@ 

$(objprefix)/double/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_DOUBLE  -c $< -o $@

$(objprefix)/complexdouble/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_DOUBLE -DLTFAT_COMPLEXTYPE -c $< -o $@

$(objprefix)/common/s%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_SINGLE -c $< -o $@

$(objprefix)/single/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_SINGLE  -c $< -o $@

$(objprefix)/complexsingle/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_SINGLE -DLTFAT_COMPLEXTYPE -c $< -o $@

$(buildprefix):
	@$(MKDIR) $(buildprefix)

$(objprefix)/common:
	@$(MKDIR) $(objprefix)$(PS)common

$(objprefix)/double:
	@$(MKDIR) $(objprefix)$(PS)double

$(objprefix)/single:
	@$(MKDIR) $(objprefix)$(PS)single

$(objprefix)/complexdouble:
	@$(MKDIR) $(objprefix)$(PS)complexdouble

$(objprefix)/complexsingle:
	@$(MKDIR) $(objprefix)$(PS)complexsingle


.PHONY: clean help doc doxy static shared munit

static: $(DTARGET) $(STARGET) $(DSTARGET)

shared: $(SO_DTARGET) $(SO_STARGET) $(SO_DSTARGET)

clean:
	@$(RMDIR) build
	@$(RMDIR) obj

help:
	@echo "USAGE: make [target]"
	@echo "Options:"
	@echo "    make [target] CONFIG=debug               Compiles the library in a debug mode"
	@echo "    make [target] NOBLASLAPACK=1             Compiles the library without BLAS and LAPACK dependencies"
	@echo "    make [target] USECPP=1                   Compiles the library using a C++ compiler"

doxy: doc
doc:
	doxygen doc/doxyconfig

cleandoc:
	@$(RMDIR) html
	@$(RMDIR) latex

munit:
	$(MAKE) clean
	# $(MAKE) BLASLAPACKLIBS="-L$(MATLABPATH) -lmwblas -lmwlapack" $(SO_DSTARGET)
	$(MAKE) $(SO_DSTARGET) COMPTARGET=debug
	$(MAKE) $(buildprefix)/ltfat.h USECPP=0

$(buildprefix)/ltfat.h: $(buildprefix) 
	$(CC) -E -P -DNOSYSTEMHEADERS -Iinclude -Ithirdparty -nostdinc include/ltfat.h -o $(buildprefix)/ltfat.h
	sed -i '1 i #ifndef _LTFAT_H' $(buildprefix)/ltfat.h
	sed -i '1 a #define _LTFAT_H' $(buildprefix)/ltfat.h
	sed -i '2 a #include <fftw3.h>' $(buildprefix)/ltfat.h
	sed -i '$$ a #endif' $(buildprefix)/ltfat.h
	$(CC) -E -P -DNOSYSTEMHEADERS -Iinclude -Ithirdparty -nostdinc $(buildprefix)/ltfat.h -o $(buildprefix)/ltfat_flat.h

install:
	install -d $(LIBDIR)
	install $(STARGET) $(DTARGET) $(DSTARGET) $(SO_STARGET) $(SO_DTARGET) $(SO_DSTARGET) $(LIBDIR)
	mkdir -p $(INCDIR)
	cp -r include/* $(INCDIR)

uninstall:
	rm -f $(LIBDIR)/$(DSTATIC) $(LIBDIR)/$(SSTATIC) $(LIBDIR)/$(DSSTATIC)
	rm -f $(LIBDIR)/$(DSHARED) $(LIBDIR)/$(SSHARED) $(LIBDIR)/$(DSSHARED)
	rm -f $(INCDIR)/ltfat.h
	rm -rf $(INCDIR)/ltfat
