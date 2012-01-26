function d=ltfatgetdefaults(fname)
%LTFATGETDEFAULTS  Get default parameters of function
%
%   `ltfatgetdefaults(fname)` returns the default parameters
%   of the function `fname` as a cell array.
%
%   `ltfatgetdefaults('all')` returns all the set defaults.
%
%   See also: ltfatsetdefaults, ltfatstart

if nargin<1
    error('%s: Too few input arguments',upper(mfilename));
end;

if strcmpi(fname,'all')
  d=ltfatarghelper('all');
else
  d=ltfatarghelper('get',fname);
end;

