function ltfatsetdefaults(fname,varargin)
%LTFATSETDEFAULTS  Set default parameters of function
%
%  LTFATSETDEFAULTS(fname,...) will set the default parameters
%  to be the parameters specified at the end of the list of input arguments.
%
%  LTFATSETDEFAULTS(fname) will clear any default parameters for the function
%  fname.
%
%  LTFATSETDEFAULTS('clearall') will clear all defaults from all
%  functions.
%
%  See also: ltfatgetdefaults, ltfatstart

if strcmpi(fname,'clearall')
  ltfatarghelper('clearall');
else
  ltfatarghelper('set',fname,varargin);
end;