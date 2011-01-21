function cout = phaseunlock(cin,a)
%PHASEUNLOCK  Undo phase lock of Gabor coefficients
%   Usage:  c=phaseunlock(c,a);
%
%   PHASEUNLOCK(c,a) removes phase locking from the Gabor coefficients c.
%   The coefficient must have been obtained from a DGT with parameter _a.
%
%   Phaselocking the coefficients modifies them so as if they were obtained
%   from a time-invariant Gabor system. A filter bank produces phase locked
%   coeffiecients. 
%
%   See also: dgt, phaselock
%
%R  puc95

%   AUTHOR:    Peter Balazs, Peter Soendergaard.
%   TESTING:   OK
%   REFERENCE: OK
  
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
phase = exp(-2*i*pi*phase);

% Handle multisignals
cout=zeros(size(cin));
for w=1:size(cin,3)
  cout(:,:,w) = cin(:,:,w).*phase;
end;

