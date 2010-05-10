function f=podd(f,dim)
%PODD   Odd part of periodic function
%   Usage:  fe=podd(f);
%           fe=podd(f,dim);
%
%   PODD(f) returns the odd part of the periodic sequence f.
%
%   PODD(f,dim) does the same along dimension dim.
%M
%   See also:  peven, dft, involute, pconv
  
f=(f-involute(f))/2;

