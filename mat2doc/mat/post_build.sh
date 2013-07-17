#!/bin/bash

cd {INST}
rm -Rf reference
rm -Rf testing
rm -Rf timing
rm -Rf mex
rm -Rf mulaclab
rm -Rf thirdparty/GPC
rm -Rf thirdparty/PolygonClip
rm mulaclab.m
mv src ..
mv thirdparty ..
mv oct ..
mv lib ..
cd ..

cp -R {CONFDIR}/src/* src/
cd src/
./bootstrap
cd ..


