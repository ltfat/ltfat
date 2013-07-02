#!/bin/bash

mat2doc mat /home/peter/nw/ltfat

cd {INST}
rm -Rf reference
rm -Rf testing
rm -Rf timing
rm -Rf mex
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
cd ..

cp -R {CONFDIR}/src .
cp {INST}/src/configure.ac src
cp {INST}/src/bootstrap src
cd ..
