function x = postpad (x, L, varargin)
%POSTPAD   Pads or truncates a vector x to a specified length L.
%   Usage: y=postpad(x,L);
%          y=postpad(x,L,C);
%          y=postpad(x,L,C,dim);
%
%   `postpad(x,L)` will add zeros to the end of the vector *x*, until the
%   result has length *L*. If *L* is less than the length of the signal, it
%   will be truncated. `postpad` works along the first non-singleton
%   dimension.
%
%   `postpad(x,L,C)` will add entries with a value of *C* instead of zeros.
%
%   `postpad(x,L,C,dim)` works along dimension *dim* instead of the first
%   non-singleton.
%
%   See also: middlepad
  
%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: NA

if nargin<2
  error('Too few input parameters.');
end

definput.keyvals.dim  = [];
definput.keyvals.C    = 0;
[~,~,C,dim] = ltfatarghelper({'C','dim'},definput,varargin,'postpad');

[x,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(x,L,dim,'POSTPAD');

if Ls<L
  x=[x; C*ones(L-Ls,W)];
else
  x=x(1:L,:);
end
  
x=assert_sigreshape_post(x,dim,permutedsize,order);
