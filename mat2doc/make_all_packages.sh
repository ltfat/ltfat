#!/bin/sh

mat2doc.py `pwd`/.. mat --script=release.py --tgz 
mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --encoding=windows-1252 --addon=ltfat_win64 --packagename=ltfat-win64
mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --addon=ltfat_maci --packagename=ltfat-maci
mat2doc.py `pwd`/.. mat --script=release.py --octpkg

