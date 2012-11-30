function xo=qam4(xi)
%QAM4  Quadrature amplitude modulation of order 4
%   Usage:  xo=qam4(xi);
%
%   `qam4(xi)` converts a vector of 0's and 1's into the complex roots of
%   unity (QAM4 modulation). Every 2 input coefficients are mapped into 1
%   output coefficient.
%
%   See also: iqam4
%
%   Demos: demo_ofdm

% Verify input

if ~all((xi==0) + (xi==1))
  error('Input vector must consist of only 0s and 1s');
end;

% Define the optimal ordering of bits
bitorder=[0;1;3;2];

% nbits is number of bits used. Everything will be ordered
% in groups of this size.
nbits=2;

L=length(xi);
symbols=L/nbits;

% nbits must divide L
if rem(symbols,1)~=0
  error('Length of input must be a multiple of 2');
end;

xi=reshape(xi,nbits,symbols);

two_power=(2.^(0:nbits-1)).';

% This could be vectorized by a repmat.
xo=zeros(symbols,1);
xo(:)=sum(bsxfun(@times,xi,two_power));

% xo now consist of numbers in the range 0:3
% Convert to the corresponding complex root of unity.
xo=exp(2*pi*i*bitorder(xo+1)/4);

