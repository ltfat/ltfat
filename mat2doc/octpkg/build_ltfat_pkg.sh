#!/bin/bash

rm -Rf ltfat
rm ltfat-1.4.0.tar.gz
mkdir ltfat
mat2doc mat /home/peter/nw/ltfat
mat2doc php /home/peter/nw/ltfat
cp -R /home/peter/publish/ltfat-mat ltfat

cd ltfat
mv ltfat-mat inst

cd inst
rm -Rf reference
rm -Rf testing
rm -Rf timing
rm -Rf mex
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
mv COPYING ..
mv CITATION ..
mv NEWS ..
cd ..

cp ~/nw/ltfat/mat2doc/octpkg/PKG_ADD .
cp ~/nw/ltfat/mat2doc/octpkg/DESCRIPTION .
cp -R ~/nw/ltfat/mat2doc/octpkg/src .
cp ~/nw/ltfat/src/configure.ac src
cp ~/nw/ltfat/src/bootstrap src
cp ~/publish/ltfat-php/INDEX .

cd ..

tar zcvf ltfat-1.4.0.tar.gz ltfat
