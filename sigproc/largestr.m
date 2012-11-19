function [xo,N]=largestr(xi,p,mtype)
%LARGESTR   Keep fixed ratio of largest coefficients
%   Usage:  xo=largestr(x,p);
%           xo=largestr(x,p,mtype);  
%           [xo,N]=largestr(...);
%
%   `largestr(x,p)` returns an array of the same size as *x* keeping
%   the fraction *p* of the coefficients. The coefficients with the largest
%   magnitude are kept.
%
%   `[xo,n]=largestr(xi,p)` additionally returns the number of coefficients
%   kept.
%
%   `largestr(x,p,'full')` returns the output as a full matrix. This is the
%   default.
%
%   `largestr(x,p,'sparse')` returns the output as a sparse matrix.
% 
%   Note that if this function is used on coefficients coming from a
%   redundant transform or from a transform where the input signal was
%   padded, the coefficient array will be larger than the original input
%   signal. Therefore, the number of coefficients kept might be higher
%   than expected.
%
%   See also:  largestn
%
%   References: ma98

%   AUTHOR : Peter L. Søndergaard
%   TESTING: OK
%   REFERENCE: OK

error(nargchk(2,3,nargin));

if (prod(size(p))~=1 || ~isnumeric(p))
  error('p must be a scalar.');
end;

if nargin==2
  mtype='full';
end;

% Determine the size of the array.
ss=prod(size(xi));

N=round(ss*p);

xo=largestn(xi,N,mtype);

