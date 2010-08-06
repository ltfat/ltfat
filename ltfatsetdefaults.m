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

global TF_CONF;

if strcmpi(fname,'clearall')
    TF_CONF.fundefs=struct;
else
    TF_CONF.fundefs.(fname)=varargin;
end;