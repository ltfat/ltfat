#This file produces system and environment dependent tools for working with files
#It sets the following variables
#
#  RM
#  CP

ifeq ($(OS),Windows_NT) 
RM = del /Q /F
CP = copy /Y
ifdef ComSpec
SHELL := $(ComSpec)
endif
ifdef COMSPEC
SHELL := $(COMSPEC)
endif
CC = gcc
else
#If not on Windows
RM = rm -rf
CP = cp -f
endif
