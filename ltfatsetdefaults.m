function ltfatsetdefaults(fname,varargin)
%LTFATSETDEFAULTS  Set default parameters of function
%
%   `ltfatsetdefaults(fname,...)` sets the default parameters to be the
%   parameters specified at the end of the list of input arguments.
%
%   `ltfatsetdefaults(fname)` clears any default parameters for the function
%   `fname`.
%
%   `ltfatsetdefaults('clearall')` clears all defaults from all functions.
%
%   See also: ltfatgetdefaults, ltfatstart

if strcmpi(fname,'clearall')
  ltfatarghelper('clearall');
else
  ltfatarghelper('set',fname,varargin);
end;

