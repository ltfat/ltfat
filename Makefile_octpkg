all:
	$(MAKE) -f Makefile_unix
	$(MAKE) -C ../oct -f Makefile_unix
	mv ../oct/*.oct .
	ls -1 *.oct | cut -d . -f 1 | xargs -i rm -f ../inst/comp/{}.m
