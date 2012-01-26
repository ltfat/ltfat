function bp = ltfatbasepath;
%LTFATBASEPATH  The base path of the LTFAT installation
%   Usage: bp = ltfatbasepath;
%
%   `ltfatbasepath` returns the top level directory in which the LTFAT
%   files are installed.
%
%   See also: ltfatstart
  
f=mfilename('fullpath');

bp = f(1:end-13);

