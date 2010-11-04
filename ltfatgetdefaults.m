function d=ltfatgetdefaults(fname)
%LTFATGETDEFAULTS  Get default parameters of function
%
%  LTFATGETDEFAULTS(fname) will return the default parameters
%  of the function fname as a cell array.
%
%  LTFATGETDEFAULTS('all') will return all the set defaults.
%
%  See also: ltfatsetdefaults, ltfatstart

if nargin<1
    error('%s: Too few input arguments',upper(mfilename));
end;

if strcmpi(fname,'all')
  d=ltfatarghelper('all');
else
  d=ltfatarghelper('get',fname);
end;