#!/bin/sh

mat2doc `pwd`/.. mat --script=release.py --tgz 
mat2doc `pwd`/.. mat --script=release.py --zip --dos --encoding=windows-1252 --addon=ltfat-win64 --packagename=ltfat-win64
mat2doc `pwd`/.. mat --script=release.py --zip --dos --addon=ltfat_maci --packagename=ltfat-maci
mat2doc `pwd`/.. mat --script=release.py --octpkg

