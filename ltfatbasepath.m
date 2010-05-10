function bp = ltfatbasepath;
%LTFATBASEPATH  The base path of the LTFAT installation
%   Usage: bp = ltfatbasepath;
%
%   LTFATBASEPATH returns the top level directory in which the LTFAT
%   files are installed.
%
%   See also: ltfatstart
  
global TF_CONF

bp = TF_CONF.basepath;