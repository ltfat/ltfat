function cout = phaseunlock(c,a)
%PHASEUNLOCK  Undo phase lock of Gabor coefficients
%   Usage:  c=phaseunlock(c,a);
%
%   `phaseunlock(c,a)` removes phase locking from the Gabor coefficients *c*.
%   The coefficient must have been obtained from a |dgt|_ with parameter *a*.
%
%   Phaselocking the coefficients modifies them so as if they were obtained
%   from a time-invariant Gabor system. A filter bank produces phase locked
%   coeffiecients. 
%
%   See also: dgt, phaselock, symphase
%
%   References: puc95

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

M=size(c,1);
N=size(c,2);
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
c=bsxfun(@times,c,phase);
