#!/bin/bash

cd {INST}
rm -Rf mex
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
mv src ..
mv oct ..
mv lib ..

# Remove Unicode characters, makeinfo in Octave cannot currenyly handle them
find -name "*.m" | xargs -n1 sed -i s/ø/oe/g
find -name "*.m" | xargs -n1 sed -i s/ö/oe/g
find -name "*.m" | xargs -n1 sed -i s/ä/a/g
find -name "*.m" | xargs -n1 sed -i s/ü/u/g
find -name "*.m" | xargs -n1 sed -i s/é/e/g
find -name "*.m" | xargs -n1 sed -i s/è/e/g
find -name "*.m" | xargs -n1 sed -i s/í/i/g

cd ..

cd src/
mv Makefile_octpkg Makefile
./bootstrap
cd ..
