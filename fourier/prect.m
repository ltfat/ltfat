function f=prect(L,n)
%PRECT   Periodic rectangle
%   Usage:  f=prect(L,n);
%
%   PSINC(L,n) computes the periodic rectangle (or square) function of
%   length L supported on n samples. The DFT of the periodic rect
%   function in the periodic sinc function.
%
%   * If n is odd, the output will be supported on exactly n samples
%   centered around the first sample.
%
%   * If n is even, the output will be supported on exactly n+1 samples
%   centered around the first sample. The function value on the two
%   samples on the edge of the function will have half the magnitude of
%   the other samples.
%
%   See also: psinc
  
error(nargchk(2,2,nargin));

if ~(numel(L)==1) || ~(isnumeric(L)) || mod(L,1)~=0 || L<=0
  error('%s: L has to be a positive integer.',upper(mfilename));
end;

if ~(numel(n)==1) || ~(isnumeric(L)) || mod(n,1)~=0 || n<=0
  error('%s: n has to be a positive integer.',upper(mfilename));
end;

f=pbspline(L,0,n);
