#!/bin/sh

make -C `pwd`/../src -f Makefile_crossmingw CROSS=/home/susnak/prg/mxe/usr/bin/x86_64-w64-mingw32- MATLIBS=/home/susnak/Dropbox/win64libs OUTDIR=/home/susnak/Dropbox/ltfat_win64 EXT=mexw64 PORTAUDIOLIB=portaudio_x64.dll
make -C `pwd`/../src -f Makefile_crossmingw CROSS=/home/susnak/prg/mxe/usr/bin/i686-pc-mingw32- MATLIBS=/home/susnak/Dropbox/win32libs OUTDIR=/home/susnak/Dropbox/ltfat_win32 EXT=mexw32 PORTAUDIOLIB=portaudio_x86.dll

mat2doc.py `pwd`/.. mat --script=release.py --tgz --unix

mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --addon=ltfat_win64 --packagename=ltfat-%s-win64
mat2doc.py `pwd`/.. mat --script=release.py --zip --dos --addon=ltfat_win32 --packagename=ltfat-%s-win32
mat2doc.py `pwd`/.. mat --script=release_keep_tests.py --octpkg --unix 

echo "Remember to push to sourceforge for the Octave package:"
echo "---> git push ssh://soender@git.code.sf.net/p/octave/ltfat master"
