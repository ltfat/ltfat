function cout = symphase(cin,a)
%SYMPHASE  Change Gabor coefficients to symmetric phase
%   Usage:  c=symphase(c,a);
%
%   SYMPHASE(c,a) alters the phase of the Gabor coefficients c so as if they
%   were obtained from a Gabor tranform based on symmetric time/frequency
%   shifts. The coefficient must have been obtained from a DGT with
%   parameter _a.
%
%   See also: dgt phaselock
%

%   AUTHORS : Peter Balazs
%             Peter Soendergaard.

error(nargchk(2,2,nargin));

if  (prod(size(a))~=1 || ~isnumeric(a))
  error('a must be a scalar');
end;

if rem(a,1)~=0
  error('a must be an integer');
end;

M=size(cin,1);
N=size(cin,2);
L=N*a;
b=L/M;

if rem(b,1)~=0
  error('Lattice error. The a parameter is probably incorrect.');
end;

TimeInd = (0:(N-1))/N;
FreqInd = (0:(M-1))*b;
    
phase = FreqInd'*TimeInd;
phase = exp(i*pi*phase);

% Handle multisignals
cout=zeros(size(cin));
for w=1:size(cin,3)
  cout(:,:,w) = cin(:,:,w).*phase;
end;
