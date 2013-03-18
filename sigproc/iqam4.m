function xo=iqam4(xi)
%IQAM4  Inverse QAM of order 4
%
%    `iqam4(xi)` demodulates a signal mapping the input coefficients to the
%    closest complex root of unity, and returning the associated bit
%    pattern. This is the inverse operation of |qam4|.
%
%    See also: qam4
%
%    Demos: demo_ofdm

% Define the optimal ordering of bits
bitorder=[0;1;3;2];

% nbits is number of bits used. Everything will be ordered
% in groups of this size.
nbits=2;

symbols=length(xi);
L=symbols*nbits;

% We round the argument of the complex numbers to the closest root of
% unity of order 4
work=round(angle(xi)/(2*pi)*4);

% work now contains negative numbers. Get rid of these
work=mod(work,4);

% Reverse the optimal bit ordering.
reversebitorder=zeros(4,1);
reversebitorder(bitorder+1)=(0:3).';

% Apply the reverse ordering
work=reversebitorder(work+1);

xo=zeros(nbits,symbols);
% Reconstruct the bits
for bit=1:nbits
  xo(bit,:)=bitget(work,bit).';
end;

xo=xo(:);
