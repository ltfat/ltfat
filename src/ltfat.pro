TEMPLATE = lib

SOURCES = dgt.c dgt_fac.c dgtreal_fac.c idgt_fac.c	\
	dgt_fb.c dgt_walnut.c gabdual_fac.c gabtight_fac.c \
	idgt_fb.c wfac.c iwfac.c \
	reassign.c windows.c \
	spread.c heapint.c winmanip.c \
	dwilt.c gabdual.c gabtight.c tfutil.c \
	pfilt.c integer_manip.c 

HEADER += ltfat.h

DEFINES += LTFAT_DOUBLE

CONFIG -= qt
CONFIG += dll

QMAKE_CFLAGS += -std=c99

INCLUDEPATH += ../thirdparty
