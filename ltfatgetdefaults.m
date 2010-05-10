function d=ltfatgetdefaults(fname)
%LTFATGETDEFAULTS  Get default parameters of function
%
%  LTFATGETDEFAULTS(fname) will return the default parameters
%  of the function fname as a cell array.
%
%  See also: ltfatsetdefaults, ltfatstart

global TF_CONF;

if isfield(TF_CONF,fname)
  d=TF_CONF.(fname);
else
  d={};
end;
