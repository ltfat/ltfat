function ltfatsetdefaults(fname,varargin)
%LTFATSETDEFAULTS  Set default parameters of function
%
%  LTFATSETDEFAULTS(fname,...) will set the default parameters
%  to be the parameters specified at the end of the list of input arguments.
%
%  LTFATSETDEFAULTS(fname) will clear any default parameters for the function
%  fname.
%
%  See also: ltfatgetdefaults, ltfatstart

global TF_CONF;

TF_CONF.(fname)=varargin;
