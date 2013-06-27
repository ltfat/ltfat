all:
	$(MAKE) -C ../inst/src -f Makefile_unix
	$(MAKE) -C ../inst/oct -f Makefile_unix
	mv ../inst/oct/*.oct .
	ls -1 *.oct | cut -d . -f 1 | xargs -i rm ../inst/comp/{}.m
