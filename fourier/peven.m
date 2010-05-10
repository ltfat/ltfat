function f=peven(f,dim)
%PEVEN   Even part of periodic function
%   Usage:  fe=peven(f);
%           fe=peven(f,dim);
%
%   PEVEN(f) returns the even part of the periodic sequence f.
%
%   PEVEN(f,dim) does the same along dimension dim.
%M
%   See also:  podd, dft, involute, pconv
  
f=(f+involute(f))/2;

