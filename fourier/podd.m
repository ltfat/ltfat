function f=podd(f,dim)
%PODD   Odd part of periodic function
%   Usage:  fe=podd(f);
%           fe=podd(f,dim);
%
%   `podd(f)` returns the odd part of the periodic sequence *f*.
%
%   `podd(f,dim)` does the same along dimension *dim*.
%
%   See also:  peven, dft, involute, pconv
 
if nargin==1 
  f=(f-involute(f))/2;
else
  f=(f-involute(f,dim))/2;
end;
