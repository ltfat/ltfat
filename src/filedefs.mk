files = dgt.o dgt_fac.o dgt_fb.o dgt_multi.o dgt_ola.o dgt_shear.o	\
  dgt_walnut.o dgtreal_fac.o dwilt.o filterbank.o heapint.o		\
  idgt_fac.o idgt_fb.o integer_manip.o iwfac.o pfilt.o reassign.o	\
  spread.o tfutil.o wfac.o windows.o winmanip.o ltfat_blaslapack.o	\
  gabdual.o gabdual_fac.o gabtight.o wavelets.o \

files_blaslapack = \
	gabdual.o gabtight.o gabdual_fac.o gabtight_fac.o

allfiles = $(files) $(files_blaslapack)

toCompile = $(files) c-safe-memalloc.o

