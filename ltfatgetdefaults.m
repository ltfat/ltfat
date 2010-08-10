function d=ltfatgetdefaults(fname)
%LTFATGETDEFAULTS  Get default parameters of function
%
%  LTFATGETDEFAULTS(fname) will return the default parameters
%  of the function fname as a cell array.
%
%  LTFATGETDEFAULTS('all') will return all the set defaults.
%
%  See also: ltfatsetdefaults, ltfatstart

global TF_CONF;

if nargin<1
    error('%s: Too few input arguments',upper(mfilename));
end;

if strcmpi(fname,'all')
    d=TF_CONF.fundefs;
else
    if isfield(TF_CONF.fundefs,fname)
        d=TF_CONF.fundefs.(fname);
    else
        d={};
    end;
end;