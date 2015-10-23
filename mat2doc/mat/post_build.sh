#!/bin/bash

cd {INST}
# We do not need any of these. They only work in Matlab
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
# Move these to the package top level
mv src ..
mv oct ..
rm -Rf mex 
# We need to keep the mex subdir. as we use some of the mex files
#mv mex ..
mv lib ..

# We moved src, we need to fix paths in some Makefiles
sed -i 's:../../src/ostools.mk:../../../src/ostools.mk:g' blockproc/java/Makefile

# Store contents from the testing and the reference directories in 
# private dir. so they do not pollute the namespace.
mkdir private
mv testing/* ./private/
mv reference/* ./private/
# Make only the test_all_ltfat.m user accessible.
mv private/test_all_ltfat.m .
# Then we do not need the directories anymore
rm -rf testing
rm -rf reference

# the thirdparty dir might contain Octave scripts as well as source code of oct and mex files
mkdir ../thirdparty
mv thirdparty/Playrec ../thirdparty
# ./thirdparty is still not empty
# rm -Rf thirdparty

# Remove Unicode characters, makeinfo in Octave cannot currently handle them
find -name "*.m" | xargs -n1 sed -i s/ø/oe/g
find -name "*.m" | xargs -n1 sed -i s/ö/oe/g
find -name "*.m" | xargs -n1 sed -i s/ä/a/g
find -name "*.m" | xargs -n1 sed -i s/ü/u/g
find -name "*.m" | xargs -n1 sed -i s/é/e/g
find -name "*.m" | xargs -n1 sed -i s/è/e/g
find -name "*.m" | xargs -n1 sed -i s/í/i/g

# Get current version
#ltfatversion=$(head -n 1 ltfat_version)
ltfatversion={VERSION}

cd ..
cd src/

# Update current version in configure.ac
sed -i -e "s/\[2.0.0\]/\[$ltfatversion\]/" configure.ac

mv Makefile_octpkg.in Makefile.in
./bootstrap
# Reported here http://savannah.gnu.org/bugs/?42278
rm -Rf autom4te.cache/
cd ..
