#!/bin/sh

mat2doc.py `pwd`/.. mat --script=release.py --tgz 

# The command below chokes on umlauts in blockana
#mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --encoding=windows-1252 --addon=ltfat_win64 --packagename=ltfat-win64
mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --addon=ltfat_win64 --packagename=ltfat-%s-win64
mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --addon=ltfat_maci --packagename=ltfat-%s-maci
mat2doc.py `pwd`/.. mat --script=release.py --octpkg
echo "Remember to push to sourceforge for the Octave package:"
echo "---> git push ssh://soender@git.code.sf.net/p/octave/ltfat master"
