function f=peven(f,dim)
%PEVEN   Even part of periodic function
%   Usage:  fe=peven(f);
%           fe=peven(f,dim);
%
%   `peven(f)` returns the even part of the periodic sequence *f*.
%
%   `peven(f,dim)` does the same along dimension *dim*.
%
%   See also:  podd, dft, involute, pconv
  
if nargin==1
  f=(f+involute(f))/2;
else
  f=(f+involute(f,dim))/2;
end;

