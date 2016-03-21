function loadlibltfat()

currdir = fileparts(mfilename('fullpath'));
libpath = [currdir, filesep, '..', filesep,'build',filesep,'libltfat.so'];
headerpath = [currdir, filesep, '..', filesep,'build',filesep,'ltfat.h'];
loadlibrary(libpath,headerpath,'mfilename','libltfatprotofile.m');





