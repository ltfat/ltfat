files = dgt.c dgtreal_fb.c dgt_multi.c dgt_ola.c dgt_shear.c	\
		dgtreal_long.c dwilt.c idwilt.c wmdct.c iwmdct.c \
		filterbank.c ifilterbank.c heapint.c heap.c wfacreal.c \
		idgtreal_long.c idgtreal_fb.c iwfacreal.c pfilt.c reassign_ti.c \
		windows.c  \
		dgt_shearola.c utils.c rtdgtreal.c \
		dgtrealwrapper.c dgtrealmp.c maxtree.c

files_complextransp = ci_utils.c ci_windows.c spread.c wavelets.c goertzel.c \
					  reassign.c gabdual_painless.c wfac.c iwfac.c \
					  dgt_long.c idgt_long.c dgt_fb.c idgt_fb.c ci_memalloc.c

files_blaslapack = ltfat_blaslapack.c gabdual_fac.c gabtight_fac.c

files_blaslapack_complextransp = gabdual.c gabtight.c

files_fftw_complextransp = dct.c dst.c

files_notypechange = memalloc.c error.c version.c argchecks.c \
					 dgtrealwrapper_typeconstant.c dgtrealmp_typeconstant.c  \
				   	 reassign_typeconstant.c wavelets_typeconstant.c \
					 integer_manip.c
