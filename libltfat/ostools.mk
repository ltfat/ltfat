#This file produces system and environment dependent tools for working with files
#It sets the following variables
#
#  RM
#  CP

ifeq ($(OS),Windows_NT) 
RM = del /Q /F
CP = copy /Y
PS2 = \\
PS = $(strip $(PS2))
ifndef SHELL
ifdef ComSpec
SHELL := $(ComSpec)
endif
ifdef COMSPEC
SHELL := $(COMSPEC)
endif
endif
CC = gcc
else
#If not on Windows
RM = rm -rf
CP = cp -f
PS = /
endif

#CC=gcc
