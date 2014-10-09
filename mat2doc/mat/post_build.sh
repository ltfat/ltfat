#!/bin/bash

cd {INST}
# We do not need any of these.
rm -Rf mex
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
# Move these to the package top level
mv src ..
mv oct ..
mv lib ..
# Store contents from the testing and the reference directories in 
# private dir. so they do not pollute the namespace.
mkdir private
mv testing/* ./private/
mv reference/* ./private/
# Make only the test_all_ltfat.m user accessible.
mv private/test_all_ltfat.m .
# Then we do not need these
rm -rf testing
rm -rf reference

mkdir ../thirdparty
mv thirdparty/Playrec ../thirdparty/Playrec
rm -Rf thirdparty

# Remove Unicode characters, makeinfo in Octave cannot currenyly handle them
find -name "*.m" | xargs -n1 sed -i s/ø/oe/g
find -name "*.m" | xargs -n1 sed -i s/ö/oe/g
find -name "*.m" | xargs -n1 sed -i s/ä/a/g
find -name "*.m" | xargs -n1 sed -i s/ü/u/g
find -name "*.m" | xargs -n1 sed -i s/é/e/g
find -name "*.m" | xargs -n1 sed -i s/è/e/g
find -name "*.m" | xargs -n1 sed -i s/í/i/g

# Get current version
ltfatversion=$(head -n 1 ltfat_version)

cd ..
cd src/

# Update current version in configure.ac
sed -i -e "s/\[2.0.0\]/\[$ltfatversion\]/" configure.ac

mv Makefile_octpkg.in Makefile.in
./bootstrap
# Reported here http://savannah.gnu.org/bugs/?42278
rm -Rf autom4te.cache/
cd ..
