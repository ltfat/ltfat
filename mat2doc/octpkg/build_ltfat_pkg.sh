#!/bin/bash

mat2doc mat /home/peter/nw/ltfat
mat2doc octpkg /home/peter/nw/ltfat

cd ~/publish/ltfat-octpkg/ltfat/inst
rm -Rf reference
rm -Rf testing
rm -Rf timing
rm -Rf mex
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
cd ..

cp -R ~/nw/ltfat/mat2doc/octpkg/src .
cp ~/nw/ltfat/src/configure.ac src
cp ~/nw/ltfat/src/bootstrap src
cp ~/publish/ltfat-octpkg/INDEX .
cp ~/publish/ltfat-octpkg/DESCRIPTION .

cd ..

tar zcvf ltfat-1.4.1.tar.gz ltfat
