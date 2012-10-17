function c = phaselock(c,a)
%PHASELOCK  Phaselock Gabor coefficients
%   Usage:  c=phaselock(c,a);
%
%   `phaselock(c,a)` phaselocks the Gabor coefficients *c*. The coefficient must
%   have been obtained from a |dgt|_ with parameter *a*.
%
%   Phaselocking the coefficients modifies them so as if they were obtained
%   from a time-invariant Gabor system. A filter bank produces phase locked
%   coefficients.
%
%   Phaselocking of Gabor coefficients correspond to the following transform:
%   Consider a signal *f* of length *L* and define $N=L/a$ and $b=L/M$.
%   The output from `c=phaselock(dgt(f,g,a,M),a)` is given by
%
%   ..             L-1 
%     c(m+1,n+1) = sum f(l+1)*exp(-2*pi*i*(m*b-n*a)*l/L)*conj(g(l-a*n+1)), 
%                  l=0  
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)e^{-2\pi il(mb-na)/L}\overline{g(l-an+1)}
%
%   where $m=0,\ldots,M-1$ and $n=0,\ldots,N-1$ and $l-an$ is computed modulo *L*.
%
%   See also: dgt, phaseunlock, symphase
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

TimeInd = (0:(N-1))*a;
FreqInd = (0:(M-1))/M;

phase = FreqInd'*TimeInd;
phase = exp(2*i*pi*phase);

% Handle multisignals
c=bsxfun(@times,c,phase);

